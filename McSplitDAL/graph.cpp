#include "graph.hpp"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>

constexpr int BITS_PER_UNSIGNED_INT (CHAR_BIT * sizeof(unsigned int));

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

dalGraph::Graph::Graph(unsigned int n) {
    this->n = n;
    label = std::vector<unsigned int>(n, 0u);
    adjmat = std::vector<std::vector<unsigned int>>(n, std::vector<unsigned int>(n, false));
    leaves = std::vector<std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<int>>>> (n);
}
void dalGraph::pack_leaves(dalGraph::Graph &g) {
    std::vector<int> deg(g.n, 0);

    for(int i = 0; i < g.n; i++)
        for(int j = 0; j < g.n; j++)
            if(i != j && g.adjmat[i][j])
                deg[i]++;

    for(int u = 0; u < g.n; u++) {
        for(int v = 0; v < g.n; v++) if(g.adjmat[u][v] && u != v && deg[v] == 1) {
            std::pair<unsigned int, unsigned int> labels(g.adjmat[u][v], g.label[v]);
            int pos = -1;
            for(int k = 0; ; k++) {
                if(k == int(g.leaves[u].size())) {
                    g.leaves[u].push_back(std::make_pair(labels, std::vector<int> ()));
                }
                if(g.leaves[u][k].first == labels) {
                    pos = k;
                    break;
                }
            }
//            assert(pos != -1);
            g.leaves[u][pos].second.push_back(v);
        }
        sort(g.leaves[u].begin(), g.leaves[u].end());
    }
}
dalGraph::Graph dalGraph::induced_subgraph(struct dalGraph::Graph& g, std::vector<int> vv) {
    dalGraph::Graph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<subg.n; j++)
            subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];//vv向量标志着顺序

    for (int i=0; i<subg.n; i++)
        subg.label[i] = g.label[vv[i]];
    return subg;
}

void add_edge(dalGraph::Graph& g, int v, int w, bool directed=false, unsigned int val=1) {
    if (v != w) {
        if (directed) {
            g.adjmat[v][w] |= val;
            g.adjmat[w][v] |= (val<<16);
        } else {
            g.adjmat[v][w] = val;
            g.adjmat[w][v] = val;
        }
    } else {
        // To indicate that a vertex has a loop, we set the most
        // significant bit of its label to 1
        g.label[v] |= (1u << (BITS_PER_UNSIGNED_INT-1));
    }
}

struct dalGraph::Graph readDimacsGraph(char* filename, bool directed, bool vertex_labelled) {
    struct dalGraph::Graph g(0);

    FILE* f;

    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    char* line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;
    int label;

    while (getline(&line, &nchar, f) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
            case 'p':
                if (sscanf(line, "p edge %d %d", &nvertices, &medges)!=2)
                    fail("Error reading a line beginning with p.\n");
                g = dalGraph::Graph(nvertices);
                break;
            case 'e':
                if (sscanf(line, "e %d %d", &v, &w)!=2)
                    fail("Error reading a line beginning with e.\n");
                add_edge(g, v-1, w-1, directed);
                edges_read++;
                break;
            case 'n':
                if (sscanf(line, "n %d %d", &v, &label)!=2)
                    fail("Error reading a line beginning with n.\n");
                if (vertex_labelled)
                    g.label[v-1] |= label;
                break;
            }
        }
    }

    if (medges>0 && edges_read != medges) fail("Unexpected number of edges.");

    fclose(f);
    return g;
}

struct dalGraph::Graph readLadGraph(char* filename, bool directed) {
    struct dalGraph::Graph g(0);
    FILE* f;

    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int w;

    if (fscanf(f, "%d", &nvertices) != 1)//dian de zongshu
        fail("Number of vertices not read correctly.\n");
    g = dalGraph::Graph(nvertices);

    for (int i=0; i<nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)//dian de du
            fail("Number of edges not read correctly.\n");
        for (int j=0; j<edge_count; j++) {
            if (fscanf(f, "%d", &w) != 1)
                fail("An edge was not read correctly.\n");
            add_edge(g, i, w, directed);
        }
    }

    fclose(f);
    return g;
}

int read_word(FILE *fp) {
    unsigned char a[2];
    if (fread(a, 1, 2, fp) != 2)
        fail("Error reading file.\n");
    return (int)a[0] | (((int)a[1]) << 8);
}

struct dalGraph::Graph readBinaryGraph(char* filename, bool directed, bool edge_labelled,
        bool vertex_labelled)
{
    struct dalGraph::Graph g(0);
    FILE* f;

    if ((f=fopen(filename, "rb"))==NULL)
        fail("Cannot open file");

    int nvertices = read_word(f);
    g = dalGraph::Graph(nvertices);

    // Labelling scheme: see
    // https://github.com/ciaranm/cp2016-max-common-connected-subgraph-paper/blob/master/code/solve_max_common_subgraph.cc
    int m = g.n * 33 / 100;
    int p = 1;
    int k1 = 0;
    int k2 = 0;
    while (p < m && k1 < 16) {
        p *= 2;
        k1 = k2;
        k2++;
    }

    for (int i=0; i<nvertices; i++) {
        int label = (read_word(f) >> (16-k1));
        if (vertex_labelled)
            g.label[i] |= label;
    }

    for (int i=0; i<nvertices; i++) {
        int len = read_word(f);
        for (int j=0; j<len; j++) {
            int target = read_word(f);
            int label = (read_word(f) >> (16-k1)) + 1;
            add_edge(g, i, target, directed, edge_labelled ? label : 1);
        }
    }
    fclose(f);
    return g;
}

struct dalGraph::Graph dalGraph::readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled) {
    struct dalGraph::Graph g(0);
    if (format=='D') g = readDimacsGraph(filename, directed, vertex_labelled);
    else if (format=='L') g = readLadGraph(filename, directed);
    else if (format=='B') g = readBinaryGraph(filename, directed, edge_labelled, vertex_labelled);
    else fail("Unknown graph format\n");
    return g;
}
