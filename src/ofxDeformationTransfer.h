
#include <string.h>
#include <ofMain.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/SparseQR>
// C++11
#include <thread>
#include <atomic>

using Eigen::DynamicSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::HouseholderQR;
using Eigen::Vector3d;
using Eigen::Vector3f;
typedef Eigen::Triplet<double> T;

namespace DT
{
	struct Triangle {
		int indices[3];
	};

	struct MeshModel {
		/* Vertices, norm vectors and triangular surfaces in the model */
		ofVec3f   *vertex;      /* array of all vertices in the model */
		ofVec3f   *normvec;     /* array of norm vectors */
		Triangle *triangle;		/* array of all triangular units */
						
		/* Number of vertices, norm vectors and triangles in the model */
		int n_vertex;
		int n_normvec;
		int n_triangle;
	};

	struct dtTransformer {
		MeshModel sourceA;   /* source reference model */
		MeshModel sourceB;   /* source deformed model */
		MeshModel targetA;   /* target reference. TBD...
								  It represents the reference model when
								  initializing the transformer object,
								  and turns into deformed target model in
								  deformation transfer solving phase. */

		map<int, int> tcdict;  /* triangle units correspondence */

		/* the deformation equation: AtA * x = c, where c = At * C */
		SparseMatrix<double> A, At, AtA;
		MatrixXd  *C, *c, x;
	};

	class ofxDeformationTransfer {

	public:
		ofxDeformationTransfer();
		~ofxDeformationTransfer();

		// load undeformed source and undeformed target (obj)
		void setReferenceModel(ofMesh *src_ref, ofMesh *trg_ref);
		Eigen::SparseLU< SparseMatrix<double> > solver;

		void transfer2TargetModel(ofMesh *sourceB, ofMesh *targetB);

	private:
		dtTransformer trans;
		MeshModel deformed_source;
		std::vector<T> vMat;
		void setDtTransferData(ofMesh *mesh, MeshModel &model);
		void setTriEdgeMatrix(MatrixXd &Va, MeshModel &model, int iTri);
		void setMatrixCol(MatrixXd &m, ofVec3f v, int iCol);
		void QRFactorize(const MatrixXd &a, MatrixXd &q, MatrixXd &r);
		void setMatrixBlock(MatrixXd &mBig, MatrixXd &mSmall, int iRow, int iCol);
		void setMatrixA(const int triIdx, std::vector<T> &v);
		void setMatrixF(const int triIdx, MatrixXd &F);
		ofVec3f calcNormal(ofVec3f v1, ofVec3f v2, ofVec3f v3);
	};
}
