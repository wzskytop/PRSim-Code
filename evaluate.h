#ifndef _EVALUATE_H_
#define _EVALUATE_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;


/*Compute Precision@k
topK1 contains the true top-k ids, computed by the Power Method or polled by the Monte Carlo Method;
realList contains the ids of the top-k highest SimRank w.r.t the query node;
*/
double calPrecision(vector<int> topK1,vector<int> realList,int k){
//topK1:gt;realList:algo
    int size = realList.size();
    int hitCount = 0;
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < topK1.size(); j++){
            if(topK1[j] == realList[i]){
		hitCount++;
                break;
            }
        }
    }
    return min(1.0, hitCount / (double) k);
}


/*
Compute AvgError@k, Precision@k for each query node in qfile.
The estimated values are in algofile.
The groundtruth or polled groundtruth is in gtfile
*/
void evaluate(string qfile,string algofile,string gtfile){
    double avg_max_error=0,avg_avg_error=0,avg_precision=0;

    vector<int> query_set;
    int query_node;
    ifstream query_file(qfile);
    while(query_file>>query_node){
       	query_set.push_back(query_node);
    }
    
    int query_num=(query_set.size()>50)?50:(query_set.size());
    int k=50;//AvgError@k, Precision@k
    for (int i = 0; i < query_num; ++i){//calculate single source simrank
        int u = query_set[i];
	cout<<i<<": "<<u<<endl;
        
	ifstream gtin(gtfile);

	vector<int> gtNodes;
        vector<double> gtvalues;
        vector<pair<double, int> > gtanswers;
        int realCnt = 0;
        int gt_tempNode;
        double gt_tempSim;    
	while(gtin>>gt_tempNode>>gt_tempSim){
            gtNodes.push_back(gt_tempNode);
            gtvalues.push_back(gt_tempSim);
            gtanswers.push_back(make_pair(gt_tempSim, gt_tempNode));
            realCnt++;
        }
        sort(gtanswers.begin(), gtanswers.end(), greater<pair<double, int> >());

	vector<int> topk_gt_Nodes;
	for(int x=0;x<k;x++){
		topk_gt_Nodes.push_back(gtanswers[x].second);
	}
	for(int y=k;y<gtanswers.size()-k;y++){
		if(gtanswers[y].first!=gtanswers[k-1].first){
			break;
		}
		else{
			topk_gt_Nodes.push_back(gtanswers[y].second);
		}
	}
	
        ifstream algoin(algofile);
        vector<int> algoNodes;
        vector<double> algovalues;
        vector<pair<double, int> > algoanswers;
        realCnt = 0;
        int algo_tempNode;
        double algo_tempSim;    
	while(algoin>>algo_tempNode>>algo_tempSim){
            algoNodes.push_back(algo_tempNode);
            algovalues.push_back(algo_tempSim);
            algoanswers.push_back(make_pair(algo_tempSim, algo_tempNode));
            realCnt++;
        }
        sort(algoanswers.begin(),algoanswers.end(),greater<pair<double, int> >());

	vector<int> topk_algo_Nodes;
        for(int x = 0; x < k; x++){
            topk_algo_Nodes.push_back(algoanswers[x].second);
        }

        double max_err = 0, avg_err = 0, precision = 0;
        int j0,jnode;
	
	for(int j=0;j<gtNodes.size();j++){
	     	double tmp_err = abs(gtanswers[j].first-algovalues[gtanswers[j].second]);
		if(max_err<tmp_err){
			max_err=tmp_err;
			j0=j;
		}
		avg_err +=tmp_err;
	}
	avg_max_error += max_err;
        avg_avg_error += avg_err / gtNodes.size();

        if(gtanswers[k-1].first > 0){
            precision = calPrecision(topk_gt_Nodes,topk_algo_Nodes,k);
        }
        else
            precision = 1;
	avg_precision += precision;
    }
    cout << "avg_avg_error="<<avg_avg_error/(double)query_num<<endl;
    cout << "avg_max_error="<< avg_max_error/(double)query_num<<endl;
    cout << "avg_precision="<<avg_precision/(double)query_num<<endl;
}

#endif
