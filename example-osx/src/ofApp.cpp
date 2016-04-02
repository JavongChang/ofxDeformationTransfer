#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){

	/* load reference model*/
	src_ref.load("horse-reference.ply");
	trg_ref.load("horse-01.ply");

	DT::ofxDeformationTransfer dTrans;

	/* set reference model */
	dTrans.setReferenceModel(&src_ref, &trg_ref);

	/* set deformed source model and transfer to target reference model */
	src_def.load("horse-reference.ply");
	trg_def.load("horse-reference.ply");
	dTrans.transfer2TargetModel(&src_def, &trg_def);
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	
	ofBackground(50);
	ofEnableDepthTest();
	
	light.enable();
	light.setPosition(cam.getPosition());

	cam.begin();

	ofPushMatrix();
	ofTranslate(-100, 100, 0);
	src_ref.draw();
	ofPopMatrix();

	ofPushMatrix();
	ofTranslate(-100, -100, 0);
	src_def.draw();
	ofPopMatrix();

	ofPushMatrix();
	ofTranslate(100, 100, 0);
	trg_ref.draw();
	ofPopMatrix();

	ofPushMatrix();
	ofTranslate(100, -100, 0);
	trg_def.draw();
	ofPopMatrix();
	cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	switch (key)
	{
		case 'f':
			ofToggleFullscreen();
			break;
		default:
			break;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
