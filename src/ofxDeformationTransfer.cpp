#include "ofxDeformationTransfer.h"

#define  EPSILON 0.0001
using namespace DT;

ofxDeformationTransfer::ofxDeformationTransfer()
{
	Eigen::initParallel();
}

ofxDeformationTransfer::~ofxDeformationTransfer()
{

}

void ofxDeformationTransfer::setReferenceModel(ofMesh *src_ref, ofMesh *trg_ref)
{
	setDtTransferData(src_ref, trans.sourceA);	/* source */
	setDtTransferData(trg_ref, trans.targetA);	/* target */

	/* make Sparse Matrix A */
	trans.A.resize(3*trans.targetA.n_triangle, trans.targetA.n_vertex);
	for (int j = 0; j < trans.targetA.n_triangle; j++) {
		setMatrixA(j, vMat);
	}

	trans.A.setFromTriplets(vMat.begin(), vMat.end());
	trans.At = trans.A.transpose();
	trans.AtA = trans.At * trans.A;
	solver.compute(trans.AtA);
	if (solver.info() != Eigen::Success){
		/* decomposition failed	*/
		cout << solver.lastErrorMessage() << endl;
	}
}

void ofxDeformationTransfer::setMatrixA(const int triIdx, std::vector<T> &v)
{
	int i1 = trans.targetA.triangle[triIdx].indices[0];
	int i2 = trans.targetA.triangle[triIdx].indices[1];
	int i3 = trans.targetA.triangle[triIdx].indices[2];
	ofVec3f e0 = trans.targetA.vertex[i2] - trans.targetA.vertex[i1];
	ofVec3f e1 = trans.targetA.vertex[i3] - trans.targetA.vertex[i1];

	/* construct Va */
	MatrixXd Va(3, 2);
	setMatrixCol(Va, e0, 0);
	setMatrixCol(Va, e1, 1);

	/* QR Decomposition */
	MatrixXd Q(3, 2); MatrixXd R(2, 2); R.setZero();
	QRFactorize(Va, Q, R);
	MatrixXd invRQT = R.inverse() * Q.transpose();

	/* set Sparse Matrix A */

	/* i1 */
	v.push_back(T(3 * triIdx + 0, i1, -invRQT(0, 0) - invRQT(1, 0)));
	v.push_back(T(3 * triIdx + 1, i1, -invRQT(0, 1) - invRQT(1, 1)));
	v.push_back(T(3 * triIdx + 2, i1, -invRQT(0, 2) - invRQT(1, 2)));

	/* i2 */
	v.push_back(T(3 * triIdx + 0, i2,  invRQT(0, 0)));
	v.push_back(T(3 * triIdx + 1, i2,  invRQT(0, 1)));
	v.push_back(T(3 * triIdx + 2, i2,  invRQT(0, 2)));

	/* i3 */
	v.push_back(T(3 * triIdx + 0, i3,  invRQT(1, 0)));
	v.push_back(T(3 * triIdx + 1, i3,  invRQT(1, 1)));
	v.push_back(T(3 * triIdx + 2, i3,  invRQT(1, 2)));
}

void ofxDeformationTransfer::transfer2TargetModel(ofMesh *sourceB, ofMesh *targetB)
{
	MatrixXd F(3 * trans.sourceA.n_triangle, 3);
	F.setZero();

	/* load deformed source model */
	setDtTransferData(sourceB, trans.sourceB);

	/* set Matrix F */
#if 1
	/* parallel computing */
	std::vector<std::thread> workers;
	std::atomic<int> i(0);

	int numThreads = std::thread::hardware_concurrency();
	numThreads = (numThreads < 1) ? 1 : numThreads;

	for (auto t = 0; t < numThreads; t++) {
		workers.push_back(std::thread([&, t]() {
			int index = 0;

			while ((index = i++) < trans.targetA.n_triangle) {
				setMatrixF(i - 1, F);	/* processing */
			}
		}));
	}

	for (auto &t : workers) {
		t.join();
	}
#else
	for (int j = 0; j < trans.sourceA.n_triangle; j++) {
		setMatrixF(j, F);
	}
#endif

	/* solve deformation transfer */
	MatrixXd UtS = trans.At * F;
	trans.x = solver.solve(UtS);

	if (solver.info() != Eigen::Success){
		/* solving failed */
		cout << solver.lastErrorMessage() << endl;
	}

	/* apply the deformation */
	for (int i = 0; i < trans.targetA.n_vertex; i++) {
		targetB->setVertex(i, ofVec3f(trans.x(i, 0), trans.x(i, 1), trans.x(i, 2)));
	}
}

void ofxDeformationTransfer::setMatrixF(const int triIdx, MatrixXd &F)
{
	/* deformation gradient for source model */
	MatrixXd Va(3, 3);
	MatrixXd Vb(3, 3);
	setTriEdgeMatrix(Va, trans.sourceA, triIdx);
	setTriEdgeMatrix(Vb, trans.sourceB, triIdx);

	/* RQ decomposition */
	MatrixXd Q(3, 3); MatrixXd R(3, 3); R.setZero();
	QRFactorize(Va, Q, R);
	MatrixXd invRQT = R.inverse() * Q.transpose();
	MatrixXd Sa = Vb * invRQT;
	Sa.transposeInPlace();
	setMatrixBlock(F, Sa, triIdx * 3, 0);
}

void ofxDeformationTransfer::setMatrixCol(MatrixXd &m, ofVec3f v, int iCol)
{
	m(0, iCol) = v.x;
	m(1, iCol) = v.y;
	m(2, iCol) = v.z;
}

void ofxDeformationTransfer::setTriEdgeMatrix(MatrixXd &Va, MeshModel &model, int iTri)
{
	int i1 = model.triangle[iTri].indices[0];
	int i2 = model.triangle[iTri].indices[1];
	int i3 = model.triangle[iTri].indices[2];

	ofVec3f v1 = model.vertex[i1];
	ofVec3f v2 = model.vertex[i2];
	ofVec3f v3 = model.vertex[i3];

	ofVec3f e1 = v2 - v1;
	ofVec3f e2 = v3 - v1;
	ofVec3f e3 = calcNormal(v1, v2, v3);

	setMatrixCol(Va, e1, 0);
	setMatrixCol(Va, e2, 1);
	setMatrixCol(Va, e3, 2);
}

void ofxDeformationTransfer::setDtTransferData(ofMesh *mesh, MeshModel &model)
{
	model.n_vertex = mesh->getNumVertices();		/* number of vertices */
	model.n_normvec = mesh->getNumNormals();		/* number of normals */
	model.n_triangle = mesh->getNumIndices() / 3;	/* number of triangles */

	model.vertex = mesh->getVerticesPointer();	/* vertex */
	model.normvec = mesh->getNormalsPointer();	/* normal */

	/* triangle */
	model.triangle = new Triangle[model.n_triangle];
	for (int i = 0; i < model.n_triangle; i++) {
		for (int j = 0; j < 3; j++) {
			model.triangle[i].indices[j] = mesh->getIndex(i * 3 + j);
		}
	}

	cout << "number of vertices : " << model.n_vertex << endl;
	cout << "number of normals : " << model.n_normvec << endl;
	cout << "number of triangles : " << model.n_triangle << endl;
	cout << "---------------------" << endl;
}

void ofxDeformationTransfer::QRFactorize(const MatrixXd &a, MatrixXd &q, MatrixXd &r)
{
	int i, j, imax, jmax;
	imax = a.rows();
	jmax = a.cols();

	for (j = 0; j<jmax; j++)
	{
		Eigen::VectorXd v(a.col(j));
		for (i = 0; i<j; i++)
		{
			Eigen::VectorXd qi(q.col(i));
			r(i, j) = qi.dot(v);
			v = v - r(i, j)*qi;
		}
		float vv = (float)v.squaredNorm();
		float vLen = sqrtf(vv);
		if (vLen < EPSILON)
		{
			r(j, j) = 1;
			q.col(j).setZero();
		}
		else
		{
			r(j, j) = vLen;
			q.col(j) = v / vLen;
		}
	}
}

ofVec3f ofxDeformationTransfer::calcNormal(ofVec3f v1, ofVec3f v2, ofVec3f v3)
{
	return v1 + ((v2 - v1).crossed(v3 - v1)) / sqrt(((v2 - v1).crossed(v3 - v1)).length());
}

void ofxDeformationTransfer::setMatrixBlock(MatrixXd &mBig, MatrixXd &mSmall, int iRow, int iCol)
{
	int r, c, rmax = mSmall.rows(), cmax = mSmall.cols();
	for (r = 0; r<rmax; r++)
	{
		for (c = 0; c<cmax; c++)
		{
			mBig(iRow + r, iCol + c) = mSmall(r, c);
		}
	}
}
