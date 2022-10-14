#include "pcatrainer.h"
#include "xmath.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
/**
*
*
*
@date:       since Aug-23-2022
@author:     Wan-Lei Zhao
@institute:  XMU, Xiamen
*
*
*
**/

using namespace std;

void help()
{
    const char *version = "Version 1.280 stable (2015-2022)";
    cout<<" Usage:\n";

    cout<<" pcatrn -pc mat -d fn (PCA analysis by SVD)\n";
    cout<<"    -pc mat\tmatrix used for training\n";
    cout<<"    -d  fn\ttrained PCA mapping matrix\n\n";

    ///cout<<"\tCodes are developped when the author is in contract with INRIA\n";
    cout<<"\tAuthor:\t\tWan-Lei Zhao, email to wlzhao@xmu.edu.cn\n";
    cout<<"\tCopyrights:\tAll rights are reserved by the atuhor\n";
    cout<<"\tVersion:\t"<<version<<endl;
}

void test()
{
    PCATrainer::test();
    ///XMath::test();

    ///cout<<sizeof(unsigned int)<<endl;

}

int main(const int argc, const char *argv[])
{
    /**
    test();
    return 0;
    /**/

    map<string, const char*> arguments;

    if(argc < 5)
    {
        help();
        return 0;
    }
    int i = 0;
    int taskType = -1;

    for(i = 1; i < argc; i += 2)
    {
        arguments.insert(pair<string, const char*>(argv[i], argv[i+1]));
        if(!strcmp(argv[i], "-pc"))
        {
           taskType = 2;
        }
    }
    if(taskType == 2)
    {
        //MissionAgent::buildPCAMat(arguments);
        PCATrainer train;
        train.pcaLearn(arguments["-pc"], arguments["-d"]);
    }else{

    }

    arguments.clear();

    return 0;
}
