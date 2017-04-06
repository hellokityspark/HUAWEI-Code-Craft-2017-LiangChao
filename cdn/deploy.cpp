#include "deploy.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include "matrix.h"
#include <map>
#include <set>
#include <stdlib.h>
#include <time.h>
using std::string;
using std::vector;
using std::size_t;

// 用来存放单条路径，及其占用的流量
struct Road {
    vector<int> nodes;
    int needs, price;
    int cost() {
        return price * needs;
    }
    Road(const vector<int>& _nodes, int _needs, int _price):
        nodes(_nodes), needs(_needs), price(_price){}
    friend std::ostream& operator<<(std::ostream& out, const Road& r) { 
        out << "price: " << r.price << " needs: " << r.needs << " ";
        for (auto n: r.nodes) out << n << " ";
        return out;
    }
};


// link 类用来记录关于消费节点的路径转移信息
struct Server {
    int position;
    vector<Road> roads;
    int needs;
    Server(int _id, int _node, int _needs): position(_node), needs(_needs){
        roads.push_back({{_id, _node}, _needs, 0});
    }
    void insert(int to, int price) {
        position = to;
        roads[0].nodes.push_back(to);
        roads[0].price += price;
    }
    void insert(vector<Road>& sub_roads, matrix<int>& Rest) {
        position = sub_roads[0].nodes.back();
        for (auto r: sub_roads) {
            assert(r.nodes.back() == position);
            for (size_t i = 0; i < r.nodes.size() - 1; ++i) {
                Rest(r.nodes[i+1], r.nodes[i]) -= r.needs;
            }
        }
        assert(std::accumulate(sub_roads.begin(), sub_roads.end(), 0, [](int o, Road r){return o + r.needs;})==needs);
        auto _copy = roads;
        roads.clear();
        for (auto& c: _copy) {
            for (auto& r: sub_roads) {
                if (r.needs >= c.needs) {
                    for (auto n: r.nodes) {
                        if (n == c.nodes.back()) continue;
                        c.nodes.push_back(n);
                    }
                    r.needs -= c.needs;
                    c.price += r.price;
                    roads.push_back(c);
                    break;
                } else {
                    auto tmp = c;
                    for (auto n: r.nodes) {
                        if (n == tmp.nodes.back()) continue;
                        tmp.nodes.push_back(n);
                    }
                    tmp.needs = r.needs;
                    tmp.price += r.price;
                    roads.push_back(tmp);
                    c.needs -= r.needs;
                    r.needs = 0;
                }
            }
        }
    }
    void merge(const Server& rhs) {
        assert(position == rhs.position);
        needs += rhs.needs;
        for (auto& r: rhs.roads) roads.push_back(r);
    }
    friend std::ostream& operator<<(std::ostream& out, const Server& server) {
        out <<"needs: "<< server.needs <<std::endl;
        for (auto& road: server.roads) out << road << std::endl;
        return out;
    }
};




// 输出
void print(vector<Server>& servers, string& s) {
    s.clear();
    s += std::to_string(std::accumulate(servers.begin(), servers.end(), 0, [](int o, const Server& s) {return o + s.roads.size();}));
    s += "\n\n";
    for (auto server: servers) {
        for (auto r: server.roads) {
            string tmp;
            for (auto iter = r.nodes.rbegin(); iter != r.nodes.rend(); ++iter) {
                tmp += std::to_string(*iter) + " ";
            }
            tmp += std::to_string(r.needs);
            tmp += "\n";
            s += tmp;
        }
    }
    while (*(s.end()-1) == '\n')
        s.erase(s.end()-1);
}

/*bool test(matrix<int>& t) {
    for (int i = 0; i < t.rows(); ++i) {
        for (int j = 0; j < t.cols(); ++j){
            if (t(i, j) < 0) return false;
        }
    }
    return true;
}*/

/*int test2(vector<Server>& servers, matrix<int>& t) {
    int total_bandwidth = 0;
    for (auto& server: servers) 
        for (auto r: server.roads)
            total_bandwidth += r.needs * (r.nodes.size() - 1);
    for (int i = 0; i < t.rows(); ++i) 
        for (int j = 0; j < t.cols(); ++j)
            total_bandwidth += t(i, j);
    std::cout<<"total bandwidth: " << total_bandwidth << std::endl;
    return total_bandwidth;
}*/

//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num, char * filename)
{
    //test(Rest());
    // 读取
    int nodes_num, links_num, servers_num;
    int cost_per_server;
    std::stringstream ss;
    // 读取第1,3行
    ss << string(topo[0]); ss >> nodes_num >> links_num >> servers_num; ss.clear();
    ss << string(topo[2]); ss >> cost_per_server; ss.clear();
    matrix<int> A(nodes_num, nodes_num);    //临接矩阵
    matrix<int> B(nodes_num, nodes_num);   //二阶临接
    vector<vector<int> > listA(nodes_num);    //A的链表版本
    matrix<double> prA(nodes_num, nodes_num);
    matrix<int> Rest(nodes_num, nodes_num); //带宽剩余
    matrix<int> Cost(nodes_num, nodes_num); //带宽费用
    //matrix<int> Allbandwidth(nodes_num, nodes_num);//总带宽
    vector<Server> servers;
    int from, to, rest, cost;
    //读取链路
    for (int i = 0; i < links_num; ++i) {
        ss << string(topo[i + 4]); ss >> from >> to >> rest >> cost; ss.clear();
        A(from, to) = 1; A(to, from) = 1;
        Rest(from, to) = rest; Rest(to, from) = rest;
        //Allbandwidth(from, to) = rest; Allbandwidth(to, from) = rest;
        Cost(from, to) = cost; Cost(to, from) = cost;
        listA[from].push_back(to); listA[to].push_back(from);
    }
    for (int i = 0; i < nodes_num; ++i) {
        for (int j = 0; j < nodes_num; ++j) {
            for (int k = 0; k < nodes_num; ++k)  
                {  
                    B(i,j)+=A(i,k)*A(k,j);  
                }  
        }
    }
    //Server一开始放在消费节点附近
    //assert(test(Rest));
    for (int i = 0; i < servers_num; ++i) {
        int _id, _node, _needs;
        ss << string(topo[i + links_num + 5]); ss >> _id >> _node >> _needs; ss.clear();
        Server server(_id, _node, _needs);
        servers.push_back(server);
    }
    //将A转化为符号pagerank的矩阵prA
    for (int j = 0; j < nodes_num; ++j) {
        double colSum = 0;
        const double pr_alpha = 0.98;// ### page_rank 参数，精度控制在maxtrix.h
        for (int i = 0; i < nodes_num; ++i) {
            if (A(i, j)) colSum += 1;
        }
        for (int i = 0; i < nodes_num; ++i) {
            if (colSum == 0) prA(i, j) = 1.0 / nodes_num;
            prA(i, j) = pr_alpha * A(i, j) / colSum + (1 - pr_alpha) / nodes_num;
        }
    }
    //计算pagerank
    vector<double> pagerank(nodes_num, double(1.0 / nodes_num));
    prA.pageRank(pagerank);
    //Step2
    const double alpha = nodes_num * 2000; //###调参
    const double N = nodes_num*1; //###调参
    vector<std::pair<int, double> > q;
    for (auto& server: servers) {
        q.clear();
        const int max_step = 40;
        int step = 0;
        while (step < max_step && pagerank[server.roads[0].nodes.back()] < 1.0 / N) {
            ++step;
            // 计算q
            int current = server.position;
            for (auto i: listA[current]) {
                double val = alpha * (pagerank[i] - pagerank[current]) - Cost(i, current) * server.needs;
                q.push_back(std::make_pair(i, val));
            }
            typedef std::pair<int, double> Pair;
            // 按照pair的second值从大到小排列
            if (q.size() > 1)
                std::sort(q.begin(), q.end(), [](const Pair& lhs, const Pair& rhs){return lhs.second > rhs.second;});
            for (auto& p: q) {
                if (Rest(p.first, current) >= server.needs) {
                    server.insert(p.first, Cost(p.first, current));
                    Rest(p.first, current) -= server.needs;
                    break;
                }
            }
        }
    }
    //assert(test(Rest));
    //合并重合服务器
    std::sort(servers.begin(), servers.end(), [](const Server& lhs, const Server& rhs){return lhs.position < rhs.position;});
    for (size_t i = 1; i < servers.size();) {
        if (servers[i - 1].position != servers[i].position) {
            ++i;
            continue;
        }
        servers[i-1].merge(servers[i]);
        servers.erase(next(servers.begin(), i));
    }
    srand(time(NULL)); 
    //Step3 合并一阶邻居
    for (int t = 0; t < 10e7; ++t) { //停止条件自己调整
        //随机选取i,j
        // assert(total_bandwidth==test2(servers, Rest));
        int ri = rand() % servers.size();
        int rj = rand() % servers.size();
        auto& si = servers[ri];
        auto& sj = servers[rj];
        int pi =  si.position;
        int pj =  sj.position;
        if (!(A(pi, pj) || (A(pi, pj) == 0 && B(pi, pj) == 1))) continue;
        //path,flow
        vector<Road> roads;
        vector<int>  choose;
        int flow = 0; // i->j能够转移的量，最多2跳
        if (A(pi, pj)) {
            roads.push_back({{pi, pj}, Rest(pj, pi), Cost(pj, pi)});
            flow += Rest(pj, pi);
        }
        for (auto k: listA[pi]) {
            if (A(k, pj) != 1) continue;
            roads.push_back({{pi, k, pj}, std::min(Rest(pj, k), Rest(k, pi)), Cost(pj, k) + Cost(k, pi)});
            flow += std::min(Rest(pj, k), Rest(k, pi));
        }
        if (flow < si.needs) continue;
        //下面将i转移到j
        //roads按照  排序
        std::sort(roads.begin(), roads.end(), [](const Road& lhs, const Road& rhs)->bool{return lhs.price < rhs.price;});
        const double beta = nodes_num * 500; // ###调参
        int cost = 0;
        int b = si.needs;
        for (size_t i = 0; i < roads.size(); ++i) {
            int real_flow = std::min(roads[i].needs, b);
            b -= real_flow;
            choose.push_back(i);
            cost += roads[i].price * real_flow;
            if (b==0) {
                roads[i].needs = real_flow;
                break;
            }
        }
        assert(std::accumulate(choose.begin(), choose.end(), 0, [&roads](int o, int c){return o + roads[c].cost();}) == cost);
        // 满足转移的条件
        if (cost - beta * (pagerank[pj] - pagerank[pi]) < cost_per_server) {
            vector<Road> chosen_roads;
            for (auto c: choose) {
                chosen_roads.push_back(roads[c]);
            }
            si.insert(chosen_roads, Rest);
            sj.merge(si);
            servers.erase(next(servers.begin(), ri));
        }
    }
    string result;
    print(servers, result);
    const char* topo_file = result.c_str();
    // 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
    write_result(topo_file, filename);
    /*int totalcost=0;
    int Allservers=servers.size();
    for (int j = 0; j < nodes_num; ++j) {
        for (int i = 0; i < nodes_num; ++i) {
           totalcost+=(Allbandwidth(j,i) - Rest(j,i))*Cost(j,i);
        }
    }
    totalcost+=Allservers*400;
    std::cout<<totalcost;*/
}

