#include "ofxDeformationTransfer.h"

#define  EPSILON 0.0001
using namespace DT;

ofxDeformationTransfer::ofxDeformationTransfer()
{

}

ofxDeformationTransfer::~ofxDeformationTransfer()
{

}

void ofxDeformationTransfer::setReferenceModel(ofMesh *src_ref, ofMesh *trg_ref)
{
	setDtTransferData(src_ref, trans.source_ref);	/* source */
	setDtTransferData(trg_ref, trans.target);		/* target */

	/* make Sparse Matrix A */
	trans.A.resize(3*trans.target.n_triangle, trans.target.n_vertex);

	for (int j = 0; j < trans.target.n_triangle; j++) {
		int i1 = trans.target.triangle[j].indices[0];
		int i2 = trans.target.triangle[j].indices[1];
		int i3 = trans.target.triangle[j].indices[2];
		ofVec3f e0 = trans.target.vertex[i2] - trans.target.vertex[i1];
		ofVec3f e1 = trans.target.vertex[i3] - trans.target.vertex[i1];

		/* construct Wj */
		MatrixXd W(3, 2);
		W.col(0) = Vector3d(e0.x, e0.y, e0.z);
		W.col(1) = Vector3d(e1.x, e1.y, e1.z);

		/* QR Decomposition */
		MatrixXd Q(3, 2); MatrixXd R(2, 2); Q.setZero(); R.setZero();
		QRFactorize(W, Q, R);
		MatrixXd invRQT = R.inverse() * Q.transpose();

		/* set Sparse Matrix A */

		/* i1 */
		trans.A.coeffRef(3*j+0, i1) = -invRQT(0, 0) - invRQT(1, 0);
		trans.A.coeffRef(3*j+1, i1) = -invRQT(0, 1) - invRQT(1, 1);
		trans.A.coeffRef(3*j+2, i1) = -invRQT(0, 2) - invRQT(1, 2);

		/* i2 */
		trans.A.coeffRef(3*j+0, i2) = invRQT(0, 0);
		trans.A.coeffRef(3*j+1, i2) = invRQT(0, 1);
		trans.A.coeffRef(3*j+2, i2) = invRQT(0, 2);

		/* i3 */
		trans.A.coeffRef(3*j+0, i3) = invRQT(1, 0);
		trans.A.coeffRef(3*j+1, i3) = invRQT(1, 1);
		trans.A.coeffRef(3*j+2, i3) = invRQT(1, 2);
	}
	trans.At = trans.A.transpose();
	trans.AtA = trans.At * trans.A;
	solver.compute(trans.AtA);
}

void ofxDeformationTransfer::transfer2TargetModel(ofMesh *deformed_src, ofMesh *deformed_trg)
{
	MatrixXd F(3 * trans.source_ref.n_triangle, 3);

	/* load deformed source model */
	MeshModel source_def;
	setDtTransferData(deformed_src, source_def);

	/* set Matrix F */
	for (int j = 0; j < trans.source_ref.n_triangle; j++) {

		/* source reference */
		int i1 = trans.source_ref.triangle[j].indices[0];
		int i2 = trans.source_ref.triangle[j].indices[1];
		int i3 = trans.source_ref.triangle[j].indices[2];
		ofVec3f e0 = trans.source_ref.vertex[i2] - trans.source_ref.vertex[i1];
		ofVec3f e1 = trans.source_ref.vertex[i3] - trans.source_ref.vertex[i1];

		/* construct V */
		MatrixXd V(3, 2);
		V.col(0) = Vector3d(e0.x, e0.y, e0.z);
		V.col(1) = Vector3d(e1.x, e1.y, e1.z);

		MatrixXd Q(3, 2); MatrixXd R(2, 2); Q.setZero(); R.setZero();
		QRFactorize(V, Q, R);
		MatrixXd invRQT = R.inverse() * Q.transpose();

		/* source deform */
		i1 = source_def.triangle[j].indices[0];
		i2 = source_def.triangle[j].indices[1];
		i3 = source_def.triangle[j].indices[2];
		e0 = source_def.vertex[i2] - source_def.vertex[i1];
		e1 = source_def.vertex[i3] - source_def.vertex[i1];

		/* construct Vd */
		MatrixXd Vd(3, 2);
		Vd.col(0) = Vector3d(e0.x, e0.y, e0.z);
		Vd.col(1) = Vector3d(e1.x, e1.y, e1.z);

		MatrixXd S = Vd * invRQT;
		SetMatrixBlock(F, S, j*3, 0);
	}

	/* apply the deformation */
	MatrixXd UtS = trans.At * F;
	trans.x = solver.solve(UtS);

	for (int i = 0; i < trans.target.n_vertex; i++) {
		deformed_trg->setVertex(i, ofVec3f(trans.x(i, 0), trans.x(i, 1), trans.x(i, 2)));
	}
}

bool ofxDeformationTransfer::loadOBJModel(string filename, MeshModel model)
{
	return true;
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
	return v1 + (v2 - v1).crossed(v3 - v1) / sqrt(((v2 - v1).crossed(v3 - v1)).length());
}

void ofxDeformationTransfer::SetMatrixBlock(MatrixXd &mBig, MatrixXd &mSmall, int iRow, int iCol)
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