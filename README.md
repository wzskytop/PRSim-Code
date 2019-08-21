# PRSim: Sublinear Time SimRank Computation on Large Power-Law Graphs
Corresponding Authors: Zhewei Wei(zhewei@ruc.edu.cn), Xiaodong He(hexiaodong_1993@ruc.edu.cn), Hanzhi Wang (wanghzccls@163.com)

Please cite our paper if you choose to use our code.
```
@inproceedings{Wei:2019:PST:3299869.3319873,
 author = {Wei, Zhewei and He, Xiaodong and Xiao, Xiaokui and Wang, Sibo and Liu, Yu and Du, Xiaoyong and Wen, Ji-Rong},
 title = {PRSim: Sublinear Time SimRank Computation on Large Power-Law Graphs},
 booktitle = {Proceedings of the 2019 International Conference on Management of Data},
 series = {SIGMOD '19},
 year = {2019},
}
```

## Tested Environment:
- Ubuntu 16.04.10
- C++ 11
- GCC 5.4.0


## Compile the code:
```
make
```


## Parameters:
- -d \<dataset\> 
- -f \<filelabel\>
- -algo \<algorithm\>
- [-e \<epsilon\> (default 0.1)]
- [-c \<damping factor\> (default 0.6)]
- [-qn \<querynum\> (default 50)]
- [-use_idx \<whether using index\> (default true,indexnum=sqrt(n))]
- [-check_dup \<check duplicate edges\> (default 0)]


## Run the example:
Before the running process, we should first generate the sorted data file in which each edge is sorted by its target node's indegree in ascending order. The sorted data file will be stored in “/dataset” named toy_graph_sorted.txt.
```
./PRSim -d dataset/toy_graph.txt -f toy_graph -algo SORT_INDEGREE_ASC
```

The running process can be divided into two steps:  
(1) Preprocessing phase: generate the index and sort them.  
(2) Querying phase: calculate the single-source SimRank.  
```
./PRSim -d dataset/toy_graph_sorted.txt -f toy_graph -algo PREPROCESS -e 0.5
./PRSim -d dataset/toy_graph_sorted.txt -f toy_graph -algo QUERY -e 0.5
```


## Instructions:
(1) datatset/: store datasets files.  
(2) query/: store query nodes files.  
(3) algo_result/: store results files.  
