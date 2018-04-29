#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 12:11:18 2018

@author: gospozha
"""

from Bio import SeqIO
from collections import defaultdict
from graphviz import Digraph 
import argparse

def parse_inputs():
     parser = argparse.ArgumentParser(description='Graph de Brujin')
     parser.add_argument('-i', '--infile' , help='Input file' , metavar='Str',
                    type=str, required=True)
     parser.add_argument('-k', '--kmersize', help = 'Kmer size', metavar = 'Int', type = int, default = 3)
     parser.add_argument('-f', '--full', help='Output = full graph', action='store_true')
     parser.add_argument('-s', '--short', help='Output = short graph', action='store_true')
     
     args = parser.parse_args()
     return args.infile, args.kmersize, args.full, args.short
    

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.n = 2
        self.coverage = 0
    
    def calc_coverage(self, c1,c2):
        self.coverage = (c1+c2)/2


class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
    def add_read(self, read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_index in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_index:(next_kmer_index+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    
    
  
#   в узлах = конкретный камер. в ребрах - камер+1, то есть ребро перекрывается с первым и вторым 
        
    def vizualize(self):
#       создаем объект - связный граф на языке DOT
        dot = Digraph()
        
#        если указан флаг вывода графа в полном формате    
        if full_flag:
#            все наши вершины будут вершинами графа
            for v in self.vertices:
                dot.node(v)
#                для каждой вершины мы создаем ребро со следующей за ней вершиной (которая есть в словаре исходящих вершин) . 
                for n in self.vertices[v].out_edges:
                    dot.edge(v,n, label = self.vertices[v].out_edges[n][0].seq)
                
#       если выбран флаг вывода графа в сокращенной форме, выводим в качестве лейблов узлов покрытие каждого камера,
#       а в качестве лейблов ребер - покрытие и длину ребра (то есть длину камера + 1)
        if short_flag:
            for v in self.vertices:
                dot.node(v, label = str(self.vertices[v].coverage))
                for n in self.vertices[v].out_edges:
                    dot.edge(v,n, label = (str(self.vertices[v].out_edges[n][0].coverage) +' '+ str(self.k+1)))
            
#       сохранение и выводод на экран     
        dot.render('test-output/graph.gv', view=True)


if __name__ == '__main__':
    in_file, size, full_flag, short_flag = parse_inputs()
    
    dataset = in_file

    k = size
    
    my_graph = Graph(k)
    

    with open(dataset, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
          
            read = str(record.seq)
            my_graph.add_read(read)
            my_graph.add_read( str(record.reverse_complement().seq) )
 
    
    for v in my_graph.vertices:
        print('Vertex: {}, coverage: {}'.format(v, my_graph.vertices[v].coverage))
        for e in my_graph.vertices[v].out_edges:
            print('out_edges: {}'.format(e))
        for e in my_graph.vertices[v].in_edges:
            print('in_edges: {}'.format(e))
            
#   вызываем функции для подсчета покрытия и визуализации       
    my_graph.calc_init_edge_coverage()
    
    my_graph.vizualize()
            
            
            
