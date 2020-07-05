#ifndef SIMSTRUCT_H
#define SIMSTRUCT_H

#define INT_MAX 32767

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>


// sort int in increasing order
bool comp1(const int &a, const int &b) {
	return a < b;
}

bool maxScoreCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

class SimStruct{
  public:
    double singlePairMulti(int u, int v, double walk_num);
    Graph g;//Class Graph
    Random R;//Class Random
    const static int NUMTHREAD = 20;
    Random Rs[NUMTHREAD];
    unsigned int seeds[NUMTHREAD];
    int vert;//the number of vertice
    string filelabel;
    double sqrtC;                                                      
    double C_value;
    double epsilon;
    double rmax;
    double nr;
    double back_nr;
    double forward_nr;
    double forward_rmax;//PRSim forward_push part's parameter 
    double backward_rmax;//PRSim backward_search's parameter
    double avg_time;
    double *H[2];//hash map,the same as U,C,UC
    int* U[2];
    int* C[1];
    int* UC[1];
    double *reserve;
    double *residue;
    double *newReserve;
    double *newResidue;
    bool* isInArray;
    int counter;
    unordered_map<int, unordered_map<int, vector<pair<int, double> > > > idx_map;//Index's data structure
    vector<int> idx_node;
    int hit_idx_count;
    int total_back_count;
    double t_idx;
    double t_back;
    SimStruct(string name, string file_lable, double eps, double c) {
        filelabel = file_lable;
        hit_idx_count = 0;
        total_back_count = 0;
        counter = 0;
        g.inputGraph(name);
        R = Random();
        vert = g.n;
        C_value = c;
        sqrtC = sqrt(C_value);
        epsilon = eps;
	rmax = (1-sqrtC)*(1-sqrtC)/25*epsilon;
	back_nr = (int)(20*log(vert)/epsilon/epsilon);
	nr = (int)(0.1*back_nr);

        forward_rmax = rmax * 2;
        backward_rmax = rmax;
        avg_time = 0;
        H[0] = new double[vert];
        H[1] = new double[vert];
        U[0] = new int[vert];
        U[1] = new int[vert];
        C[0] = new int[vert];
        UC[0] = new int[vert];
        reserve = new double[vert];
        residue = new double[vert];
        newReserve = new double[vert];
        newResidue = new double[vert];
        isInArray = new bool[vert];
        for(int i = 0; i < vert; i++){
            isInArray[i] = false;
            H[0][i] = 0;
            H[1][i] = 0;
            U[0][i] = 0;
            U[1][i] = 0;
            C[0][i] = 0;
            UC[0][i] = -1;
        }
        srand(unsigned(time(0)));
        for(int i  =0; i < NUMTHREAD; i++){
            seeds[i] = unsigned(rand());
            Rs[i] = Random(seeds[i]);
        }
        t_idx = 0;
        t_back = 0;
        cout << "====init done!====" << endl;
    }
    ~SimStruct() {
        delete[] H[0];
        delete[] H[1];
        delete[] U[0];
        delete[] U[1];
        delete[] C[0];
        delete[] UC[0];
        delete[] reserve;
        delete[] residue;
        delete[] newReserve;
        delete[] newResidue;
        delete[] isInArray;
    }

      
    //Preprocessing: load index
    void loadIndex(string idxname){
        cout<<"idxname:"<<idxname<<endl;
	ifstream index_file;
	index_file.open(idxname.data());
	assert(index_file.is_open());
	cout << "===start load index===" << endl;
        
	clock_t t0 = clock();
        int indexsize=0;
	int u, w, level;
        double reserve;   
	while(index_file>>u>>w>>level>>reserve){
	    indexsize+=1;
            idx_map[w][level].push_back(pair<int, double>(u, reserve));
        }
        index_file.close();
	cout<<"indexsize="<<indexsize<<endl;
        cout << "====load index done!====" << endl;
        clock_t t1 = clock();
        cout << "load index time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
    } 

    bool isInIndex(int w, int l){
        if(idx_map.find(w) != idx_map.end()){
            if(idx_map[w].find(l) != idx_map[w].end()){
                for(int i = 0; i < idx_map[w][l].size(); i++){
                    cout << idx_map[w][l][i].first << " " << idx_map[w][l][i].second << endl;
                }
                return true;
            }
        }
        return false;
    }

    //Estimate pi(u,w)
    unordered_map<int, vector<pair<int, double> > > foraMap(int u){
        unordered_map<int, unordered_map<int, double> > answer;
        unordered_map<int, vector<pair<int, double> > > realAnswer;

        residue[u] = 1;
        vector<int> candidate_set;
        candidate_set.push_back(u);
        isInArray[u] = true;
        vector<int> new_candidate_set;
        int tempLevel = 0;

        clock_t t0 = clock();
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;        
        while(candidate_set.size() > 0){
            for(int j = 0; j < candidate_set.size(); j++){
                isInArray[candidate_set[j]] = false;
            }
            int residue_count = 0;
            for(int j = 0; j < candidate_set.size(); j++){
                int tempNode = candidate_set[j];
                double tempR = residue[tempNode];
                newReserve[tempNode] += (1 - sqrtC) * tempR;
                int inSize = g.getInSize(tempNode);
                for(int k = 0; k < inSize; k++){
                    int newNode = g.getInVert(tempNode, k);
                    newResidue[newNode] += residue[tempNode] * sqrtC / (double)inSize;
                    if(U[1][newNode] == 0){
                        U[0][residue_count++] = newNode;
                        U[1][newNode] = 1;
                    }
                    if(!isInArray[newNode] && newResidue[newNode] > forward_rmax){
                        isInArray[newNode] = true;
                        new_candidate_set.push_back(newNode);
                    }
                }
                residue[tempNode] = 0;  
            }
            for(int j = 0; j < candidate_set.size(); j++){
                int tempNode = candidate_set[j];
                if(newReserve[tempNode] > 0){
                    answer[tempLevel][tempNode] = newReserve[tempNode];
                }
                newReserve[tempNode] = 0;
            }
            for(int j = 0; j < residue_count; j++){
                if(newResidue[U[0][j]] <= forward_rmax){
                    rsum += newResidue[U[0][j]];
                    aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel + 1,U[0][j]), newResidue[U[0][j]]));
                }
                else{
                    residue[U[0][j]] = newResidue[U[0][j]];
                }
                newResidue[U[0][j]] = 0;
                U[1][U[0][j]] = 0;
            }
            candidate_set = new_candidate_set;
            new_candidate_set.clear();
            tempLevel++;
        }
        Alias alias = Alias(aliasP);
        clock_t t1 = clock();
        cout << "fora" << endl;
        if(rsum > 0){
            int fora_walk_num = (t1 - t0) / (double) CLOCKS_PER_SEC * 400000;
	    cout<<"fora_walk_num="<<fora_walk_num<<endl;	    
	    double increment = rsum * (1 - sqrtC) / (double)fora_walk_num;
            for(int j = 0; j < fora_walk_num; j++){
                pair<int, int> tempPair = alias.generateRandom(R);
                int tempLevel = tempPair.first;
                int tempNode = tempPair.second;
                int tempCount = 0;
                if(answer.find(tempLevel) != answer.end() && answer[tempLevel].find(tempNode) != answer[tempLevel].end())
                    answer[tempLevel][tempNode] += increment;
                else
                    answer[tempLevel][tempNode] = increment;
                while(R.drand() < sqrtC){
                    int length = g.getInSize(tempNode);
                    if(length > 0){   
                        int r = R.generateRandom() % length;
                        tempNode = g.getInVert(tempNode, r);
                        tempCount++;
                    }
                    else{
                        break;
                    }
                    if(answer.find(tempLevel + tempCount) != answer.end() && answer[tempLevel + tempCount].find(tempNode) != answer[tempLevel + tempCount].end())
                        answer[tempLevel + tempCount][tempNode] += increment;
                    else
                        answer[tempLevel + tempCount][tempNode] = increment;
                } 
            }
        }             
        clock_t t2 = clock();
        for(auto& c: answer){
            for(auto& d: c.second){
                if(d.second > forward_rmax){
                    realAnswer[c.first].push_back(pair<int, double>(d.first, d.second));
                }
            }
        } 
        return realAnswer;
    }

    
    int sampleD(int nodeId, int walk_num){
        int tempCount = 0;
        for(int i = 0; i < walk_num; i++){
            int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
            double meet_level = 0;
            while(R.drand() < C_value){
                int length = g.getInSize(u_newNode);
                if(length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if(length == 0)
                    break;
		r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);
                meet_level++;
                if(u_nextNode == v_nextNode){
                    if(meet_level > 1)
                        tempCount += 1;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        return tempCount;
    }
    
        
    //Calculating pi(v,w) using sequential probe method
    void getRandomBackList(int w, int targetLevel, int backWalkNum, double forwardSim, double dw, double* answer){
        total_back_count++;
        unordered_map<int, double> walkMap;
        bool isIdx = false;
        clock_t t0 = clock();
        if(idx_map.find(w) != idx_map.end() && idx_map[w].find(targetLevel) != idx_map[w].end()){
            isIdx = true;
	    hit_idx_count++;
            for(auto &c: idx_map[w][targetLevel]){
                if(c.second > backward_rmax){
                    answer[c.first] += dw * forwardSim * c.second / (1 - sqrtC) / (1 - sqrtC);
                }
                else
                    break;
            }
        }
        clock_t t1 = clock();
        int edge_count = 0;
        if((!isIdx)&&(find(idx_node.begin(),idx_node.end(),w)==idx_node.end())){
	    double increment = 1 / (double) backWalkNum * pow(sqrtC, targetLevel) * (1 - sqrtC);
            for(int j = 0; j < backWalkNum; j++){
                int ind = 0;
                H[ind][w] = 1;
                int Ucount = 1;
                int Ucount1 = 0;
                int UCcount = 0;
                U[0][0] = w;
                for (int i = 0; i < targetLevel; i++){
                    for (int j = 0; j < Ucount; j++){
                        int tempNode = U[ind][j];
                        int outCount = g.getOutSize(tempNode);
                        double tempMaxInSize = 1 / R.drand();
                        for (int k = 0; k < outCount; k++){
                            int newNode = g.getOutVert(tempNode, k);
                            if(g.getInSize(newNode) > tempMaxInSize){
				break;
                            }
                            if(H[1-ind][newNode] == 0){
                                U[1-ind][Ucount1++] = newNode;
                            }
                            H[1-ind][newNode] += H[ind][tempNode];
                        }
                    }
                    for (int j = 0; j < Ucount; j++){
                        H[ind][U[ind][j]] = 0;
                        U[ind][j] = 0;
                    }
                    Ucount = Ucount1;
                    Ucount1 = 0;
                    ind = 1 - ind;
                    if (Ucount == 0)
                        break;
                }
                for (int i = 0; i < Ucount; i++){
                    int tempNode = U[ind][i];
                    if(walkMap.find(tempNode) == walkMap.end()){
                        walkMap[tempNode] = H[ind][tempNode] * increment;
                    }
                    else{
                        walkMap[tempNode] += H[ind][tempNode] * increment;
                    }
                    U[ind][i] = 0;
                    H[ind][tempNode] = 0;
                }
                Ucount = 0;

            }
            for(auto &key: walkMap){
                int tempNode = key.first;
                double tempSim = key.second;
                if(tempSim > backward_rmax){
                    answer[tempNode] += dw * forwardSim * tempSim / (1 - sqrtC) / (1 - sqrtC);
                }
            }
        }
        
        clock_t t2 = clock();
        t_idx += t1 - t0;
        t_back += t2 - t1;
    }


    void generateSortByIndegreeAsc(string outputFile){
        ofstream fout(outputFile);
	fout << vert << "\n";
        for(int i = 0; i < vert; i++){
	    vector<pair<int, double> > tempVec;
            for(int j = 0; j < g.getOutSize(i); j++){
                tempVec.push_back(pair<int, double>(g.getOutVert(i, j), g.getInSize(g.getOutVert(i, j))));
            }
            sort(tempVec.begin(), tempVec.end(), maxScoreCmp);
            for(int j = tempVec.size() - 1; j >= 0; j--){
		if(i<1000){
			cout << i <<" node="<<tempVec[j].first<<" insize="<<tempVec[j].second<<endl;
		}
		fout << i << " " << tempVec[j].first << "\n";
            }
        }
        fout.close();
    }

    
    void singleSourceQuery(int u){
	cout<<"u="<<u<<endl;
	double single_walk_num=nr;
	cout<<"sampleD_walk_num="<<single_walk_num<<endl;
	cout<<"backward_walk_num="<<back_nr<<endl;

        double* answer = new double[vert];
        for(int i = 0; i < vert; i++){
            answer[i] = 0;
            isInArray[i] = false;
        }
        double forward_time = 0, alias_time=0,backward_time = 0, calD_time = 0, calAns_time = 0;
        cout << "init done!" << endl;
        clock_t t0 = clock();
        unordered_map<int, vector<pair<int, double> > > forward_map = foraMap(u);
        clock_t t1 = clock();
        forward_time = t1 - t0;

        clock_t t_alias_s = clock();
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        for(auto& c: forward_map){
            int tempLevel = c.first;
            for(auto& d: c.second){
                int w = d.first;
                double forwardSim = d.second;
                rsum += forwardSim;
                aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel,w), forwardSim));
            }
        }
        Alias alias = Alias(aliasP);
        clock_t t_alias_e = clock();
	alias_time=t_alias_e-t_alias_s;
        cout << "create Alias done! rsum: " << rsum << endl;

	unordered_map<int, int> hitMap;
        unordered_map<int, int> totalMap;

	for(int i = 0; i<(int)single_walk_num;i++){
            pair<int, int> tempPair = alias.generateRandom(R);
	    int tempLevel = tempPair.first;
            int tempNode = tempPair.second;
            int hitCount = sampleD(tempNode,1); 
	    //cout<<"hitCount="<<hitCount<<endl;
	    if(hitMap.find(tempNode) == hitMap.end())
                hitMap[tempNode] = 0;
            if(totalMap.find(tempNode) == totalMap.end())
                totalMap[tempNode] = 0;
            hitMap[tempNode] += hitCount;
	    totalMap[tempNode] +=1;
        }

        unordered_map<int, double> dValue;
        int cnt_w=0;
	for(auto& w : hitMap){
            cnt_w+=1;
	    if(g.getInSize(w.first) == 0)
                dValue[w.first] = 1;
            else 
                dValue[w.first] = 1 - C_value / (double) g.getInSize(w.first) - w.second / (double)totalMap[w.first];
        }
        clock_t t_calD_e = clock();
        cout << "cal D done!" << endl;
        calD_time += t_calD_e - t_alias_e;

	for(auto& c: forward_map){
            int tempLevel = c.first;
            for(auto& d: c.second){
                int w = d.first;
                double forwardSim = d.second;
                double dw = 0;
                if(g.getInSize(w) == 0)
                    dw = 1;
                else if(dValue.find(w) != dValue.end())
                    dw = dValue[w];
                else
                    dw = 1 - C_value / (double) g.getInSize(w);
                //cout<<"d["<<w<<"]="<<dw<<endl;
		clock_t t2 = clock();
                getRandomBackList(w, tempLevel, back_nr * forwardSim / rsum + 1, forwardSim, dw, answer); 
                clock_t t3 = clock();
                
                backward_time += t3 - t2;
            }
        }

        clock_t t_end = clock();
        cout << "cal Answer done!" << endl;
     
	answer[u]=1.0;
   
	stringstream ssout;
	ssout<<"algo_result/"<<filelabel<<"_"<<u<<"_"<<epsilon<<".txt";
	ofstream fout(ssout.str());
        vector<pair<int, double>> pprs;
        for(int j = 0; j < vert; j++){
            pprs.push_back(pair<int, double>(j, answer[j]));
	    fout<<j<<" "<<answer[j]<<"\n";
        }

        avg_time += (t_end - t0) / (double) CLOCKS_PER_SEC;
        delete[] answer;   
        return;
    }

   
    //Generate index
    void backwardSearch(int w, ofstream &idx_file, double index_rmax){
        residue[w] = 1;
        double total_r = 1, total_pi = 0;
        vector<int> candidate_set;
        candidate_set.push_back(w);
        isInArray[w] = true;
        vector<int> new_candidate_set;
        int tempLevel = 0;
        double edgeNum = 0;
        
        while(candidate_set.size() > 0){
            for(int j = 0; j < candidate_set.size(); j++){
                isInArray[candidate_set[j]] = false;
            }
            int residue_count = 0;
            for(int j = 0; j < candidate_set.size(); j++){
                int tempNode = candidate_set[j];
                double tempR = residue[tempNode];
                newReserve[tempNode] += (1 - sqrtC) * tempR;        
                int outSize = g.getOutSize(tempNode);
                double tempEdgeCount = 0;
                for(int k = 0; k < outSize; k++){
                    int newNode = g.getOutVert(tempNode, k);
                    double incre = residue[tempNode] * sqrtC / (double)g.getInSize(newNode);
                    tempEdgeCount++;
                    newResidue[newNode] += incre;
                    if(U[1][newNode] == 0){
                        U[0][residue_count++] = newNode;
                        U[1][newNode] = 1;
                    }
                    if(!isInArray[newNode] && newResidue[newNode] > index_rmax){
                        isInArray[newNode] = true;
                        new_candidate_set.push_back(newNode);
                    }
                }   
                edgeNum += tempEdgeCount;
                residue[tempNode] = 0;
                total_r--;
            }
            for(int j = 0; j < candidate_set.size(); j++){
                int tempNode = candidate_set[j];
                if(newReserve[tempNode] > index_rmax){
                    total_pi++;
                    idx_map[w][tempLevel].push_back(pair<int,double>(tempNode,newReserve[tempNode]));
		}
	        newReserve[tempNode] = 0;
            }
            for(int j = 0; j < residue_count; j++){
                if(newResidue[U[0][j]] > index_rmax)
                    residue[U[0][j]] = newResidue[U[0][j]];
                newResidue[U[0][j]] = 0;
                U[1][U[0][j]] = 0;
            }
       
	    candidate_set.swap(new_candidate_set);
	    vector<int>().swap(new_candidate_set);
          
	    sort(idx_map[w][tempLevel].begin(),idx_map[w][tempLevel].end(),maxScoreCmp);

	    tempLevel++;
        }
	
    }


    //Calculate PageRank
    void PageRankPowerMethod(int iterations, string outputFile){
        double pagerank_eps = 0.000000001;
        double* scores = new double[vert];
        double* nextScores = new double[vert];
        double initValue = 1 / (double) vert;
        double alpha = 1 - sqrtC;
        for(int i = 0; i < vert; i++){
            scores[i] = initValue;
            nextScores[i] = 0;
            
        }
        for(int i = 1; i <= iterations; i++){
            bool isConverge = true;
            for(int j = 0; j < vert; j++){
                nextScores[j] = alpha * initValue;
                for(int k = 0; k < g.getOutSize(j); k++){
                    int tempNode = g.getOutVert(j, k);
                    nextScores[j] += (1-alpha) * scores[tempNode] / (double) g.getInSize(tempNode);
                }
            }
            for(int j = 0; j < vert; j++){
                if(fabs(scores[j] - nextScores[j]) > pagerank_eps){
                    isConverge = false;
                }
                scores[j] = nextScores[j];
            }
            if(isConverge)
                break;
        }
        ofstream fout(outputFile);
        vector<pair<int, double> > pprs;
        for(int j = 0; j < vert; j++){
            pprs.push_back(pair<int, double>(j, scores[j]));
        }
        sort(pprs.begin(), pprs.end(), maxScoreCmp);
        for(int j = 0; j < vert; j++){
            if(pprs[j].second >= 0){
                fout << pprs[j].first << " " << pprs[j].second << "\n";
            }
        }
        fout.close();
        delete[] scores;
        delete[] nextScores;
    }
       
};


#endif
