#include <opencv.hpp>
#include <dirent.h>
#include <iostream>

using namespace cv;

void displayDepth(const string& name, Mat& depthMap) {
    Mat show;
    const float scaleFactor = 0.05f;
    depthMap.convertTo(show, CV_8UC1, scaleFactor);
    imshow(name, show);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "You need to supply path for orginal image folder and the recreated image folder." << std::endl;
    }

	Mat cvMat1 = Mat::zeros(480, 640, CV_16UC1);
	Mat cvMat2 = Mat::zeros(480, 640, CV_16UC1);
	std::string dirname1(argv[1]);
	std::string dirname2(argv[2]);

    DIR* dir = opendir(dirname1.c_str());
    if (dir == NULL) {
        std::cout << "Can not read directory." << std::endl;
        return -1;
    }
    std::vector<std::string> entries1;
    struct dirent* ent;
    while ((ent = readdir(dir)) != NULL) {
        entries1.push_back(std::string(ent->d_name));
    }
    closedir(dir);

    dir = opendir(dirname2.c_str());
    if (dir == NULL) {
        std::cout << "Can not read directory." << std::endl;
        return -1;
    }
    std::vector<std::string> entries2;
    while ((ent = readdir(dir)) != NULL) {
        entries2.push_back(std::string(ent->d_name));
    }
    closedir(dir);

    std::sort(entries1.begin(), entries1.end());
    std::sort(entries2.begin(), entries2.end());

	std::string depthName("depth");
    for (int i = 2; i < entries2.size(); ++i) {
    	std::string file1 = entries1[i];
    	std::string file2 = entries2[i];
    	std::cout << dirname1 + '/' + file1 << std::endl;
    	std::cout << dirname2 + '/' + file2 << std::endl;
        cvMat1 = imread(dirname1 + '/' + file1, -1);
        cvMat2 = imread(dirname2 + '/' + file2, -1);
        imshow("Original", 16*cvMat1);
        Mat diff = 100*(cvMat2 - cvMat1);
        displayDepth("First", diff);
        //displayDepth("Second", cvMat2);
        waitKey();
    }
}
