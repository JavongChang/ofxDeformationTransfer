
#include <string.h>
#include <ofMain.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/SparseQR>

using Eigen::DynamicSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::HouseholderQR;
using Eigen::Vector3d;
using Eigen::Vector3f;

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
		MeshModel source_ref;   /* source reference model */
		MeshModel target;       /* target reference/deformed model.
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

		void transfer2TargetModel(ofMesh *deformed_src, ofMesh *deformed_trg);

	private:
		dtTransformer trans;
		MeshModel deformed_source;

		// set triangle j with indices(i1, i2, i3);
		bool loadOBJModel(string filename, MeshModel model);
		void setDtTransferData(ofMesh *mesh, MeshModel &model);
		void QRFactorize(const MatrixXd &a, MatrixXd &q, MatrixXd &r);
		ofVec3f calcNormal(ofVec3f v1, ofVec3f v2, ofVec3f v3);
		void SetMatrixBlock(MatrixXd &mBig, MatrixXd &mSmall, int iRow, int iCol);
	};
}