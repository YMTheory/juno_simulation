/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Tue Sep 15 11:43:19 2020
 @File Name: ap_toy.cc
 ************************************************************************/

void ap_toy()
{
    double dynComp1_prob   = 0.591/100*6.57;
    double dynComp2_prob   = 0.276/100*1.99;
    double dynComp1_tmu    = 4.03e3;
    double dynComp1_tsigma = 1.48e3;
    double dynComp2_tmu    = 15.54e3;
    double dynComp2_tsigma = 30.1e3;

    double mcpComp1_prob   = 0.152/100*15.5;
    double mcpComp2_prob   = 0.145/100*11.4;
    double mcpComp3_prob   = 0.135/100*9.1;
    double mcpComp1_tmu    = 4.03e3;
    double mcpComp1_tsigma = 1.48e3;
    double mcpComp2_tmu    = 15.54e3;
    double mcpComp2_tsigma = 30.1e3;


    double sample = 0;
    TH1D* h1 = new TH1D("h1", "", 10, 0, 10);
    for(int i=0; i<1000000; i++) {
        sample = gRandom->Poisson(dynComp1_prob);
        cout << sample << endl;
    }

    h1->Draw();

}

