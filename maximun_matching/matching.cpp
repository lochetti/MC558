/*
	Implementation of the Hungarian Algorithm to find a maximum matching on a Km,n graph
    Copyright (C) 2012  Renato Tadeu Lochetti

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/matching.h>
#include <lemon/adaptors.h>
#include <fstream>
#include <string>

void readGraph(std::string filename, lemon::ListGraph& g,
				lemon::ListGraph::EdgeMap<int>& weight, 
				lemon::ListGraph::NodeMap<int>& partition){

	// Read the file.
	std::ifstream file(filename.c_str());
	std::stringstream input;
	if(file){
		input << file.rdbuf();
		file.close();
	} else {
		std::cout << "File not found: " << filename << std::endl;
		exit(1);
	}
	
	// Parser.
	lemon::ListGraph::NodeMap<std::string> rotulo(g);
	lemon::GraphReader<lemon::ListGraph>(g, input).
		edgeMap("weight", weight).
		nodeMap("partition", partition).
		run();
}

lemon::ListGraph::Node getNodeFromMap(lemon::ListGraph& g, lemon::ListGraph::NodeMap<int>& map, int value, lemon::ListGraph::NodeMap<int>& partition, int parte){
	for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n){
		if( map[n] == value && partition[n] == parte){
			return n;
		}
	}
	return lemon::INVALID;
}

int main(int argc, char *argv[]){
	
	if(argc != 2){
		std::cout << "uso: " << argv[0] << " <filename>" << std::endl;
		return 1;
	}

	// Read the graph.
	std::string filename = argv[1];	
	
	lemon::ListGraph g;
	lemon::ListGraph::NodeMap<int> nodeLabelMenor(g);
	lemon::ListGraph::EdgeMap<int> weight(g);
	lemon::ListGraph::NodeMap<int> partition(g);
	lemon::ListGraph::NodeMap<int> cobertura(g);
	lemon::ListGraph::NodeMap<int> nodeLabelMaior(g);
	int countParticao1 = 0, countParticao2 = 0;
	readGraph(filename, g, weight, partition);
	
	for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n){
		if(partition[n] == 1){ //one side of the partition
			countParticao1++;
		}else{ //another side
			countParticao2++;
		}
	}
	int menorParticao = countParticao1 < countParticao2 ? 1 : 0;
	int maiorParticao = menorParticao == 1 ? 0 : 1;
	int countMenorParticao = countParticao1 < countParticao2 ? countParticao1 : countParticao2;
	int countMaiorParticao = countMenorParticao == countParticao1 ? countParticao2 : countParticao1;
	int labelMenorCount = 0, labelMaiorCount = 0;
	//generate the initial matching putting the edge weight in the nodes of the small part.
	for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n){
		if(partition[n] == menorParticao){
			int peso = -9999999;
			for(lemon::ListGraph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e){
				int ePeso = weight[e];
				if(ePeso > peso)
					peso = ePeso;
			}
			cobertura[n] = peso;
			nodeLabelMenor[n] = labelMenorCount;
			
			labelMenorCount++;
		}else{
			cobertura[n] = 0;
			nodeLabelMaior[n] = labelMaiorCount;
			labelMaiorCount++;
		}
	}
	
	int **matrizDeExcesso = new int*[countMenorParticao];
	for(int i = 0; i < countMenorParticao; i++)
		matrizDeExcesso[i] = new int[countMaiorParticao];
		
	int pesoDoEmparelhamentoMaximo = 0;
	
	//main looping . Stay here while the maximum matching was not found.
	while(true){
		std::vector<lemon::ListGraph::Node> coberturaLinha;
		coberturaLinha.erase(coberturaLinha.begin(), coberturaLinha.end());
		std::vector<lemon::ListGraph::Node> r;
		r.erase(r.begin(), r.end());
		std::vector<lemon::ListGraph::Node> t;
		t.erase(t.begin(), t.end());
		lemon::ListGraph::NodeMap<bool> node_filter(g, true);
		lemon::ListGraph::EdgeMap<bool> edge_filter(g, true);
		lemon::SubGraph<lemon::ListGraph> subGrafoDeIgualdade(g, node_filter, edge_filter);
		
		//Calculation of the excess matrix and generation of the equality subgraph.
		for(int i = 0; i < countMenorParticao; i++){
			lemon::ListGraph::Node noParticaoMenor = getNodeFromMap(g, nodeLabelMenor, i, partition, menorParticao);
			if(noParticaoMenor != lemon::INVALID){
				for(int j = 0; j < countMaiorParticao; j++){
					lemon::ListGraph::Node noParticaoMaior = getNodeFromMap(g, nodeLabelMaior, j, partition, maiorParticao);
					if(noParticaoMaior != lemon::INVALID){
						lemon::ListGraph::Edge aresta = findEdge(g, noParticaoMenor, noParticaoMaior);
						int valorMatriz = cobertura[noParticaoMenor] + cobertura[noParticaoMaior] - weight[aresta];
						matrizDeExcesso[i][j] = valorMatriz;
						if(valorMatriz != 0){
							subGrafoDeIgualdade.disable(aresta);
						}
					}
				}
			}
		}
		//caltulation of the maximum matching in the equality subgraph.
		lemon::MaxMatching< lemon::SubGraph<lemon::ListGraph> > subMatching(subGrafoDeIgualdade);
		subMatching.run();
		int matchingSize = subMatching.matchingSize();
		if(matchingSize == countMenorParticao){ //found one perfect matching
			for(lemon::ListGraph::EdgeIt e(g); e!=lemon::INVALID; ++e){
				if(subMatching.matching(e)){
					pesoDoEmparelhamentoMaximo+= weight[e];
				}
			}
			break;
		}else{ //perfect matching not found . We need to re-calculate the the covers,the excess matrix and try again!
			int aux = 0;
			//cover calculation on the equality subgraph.
			for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID && matchingSize != aux ; ++n, aux++){
				if( partition[n] == maiorParticao){
					lemon::ListGraph::Edge arestaNoEmparelhamento = subMatching.matching(n);
					if(arestaNoEmparelhamento != lemon::INVALID){
						coberturaLinha.push_back(n);
					}
				}
			}
			//calculation of the node sets R and T that will be used to recalculate the vertex cover on the initial graph.
			for(std::vector<lemon::ListGraph::Node>::size_type i = 0; i != coberturaLinha.size(); i++){
				lemon::ListGraph::Node myNode = coberturaLinha[i];
				if(partition[myNode] == menorParticao){
					//add on R
					r.push_back(myNode);
				}else{
					//add on T
					t.push_back(myNode);
				}
			}
			//calculate epsilon to remake the vertex cover on the initial graph.
			int epsilon = 999999999; //something really big here.
			for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n){
				if(partition[n] == menorParticao && std::find(r.begin(), r.end(), n) == r.end()){ // n belongs to the small part and is not in R
					for(lemon::ListGraph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e){
						lemon::ListGraph::Node nVizinho = g.oppositeNode(n, e);
						if(partition[nVizinho] == maiorParticao && std::find(t.begin(), t.end(), nVizinho) == t.end()){ // nVizinho belongs to the bigger part and is not in T .
							if(weight[e] < epsilon)
								epsilon = weight[e];
						}
					}
				}
			}
			// remake the vertex cover on the initial graph
			for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n){
				if(partition[n] == menorParticao && std::find(r.begin(), r.end(), n) == r.end()){
					cobertura[n] -= epsilon;
				}
			}
			for(std::vector<lemon::ListGraph::Node>::size_type i = 0; i != t.size(); i++) {
				lemon::ListGraph::Node n = t[i];
				cobertura[n] += epsilon;
			}
		}
	}
	
	std::cout << "matchingWeight " << pesoDoEmparelhamentoMaximo << std::endl;
}