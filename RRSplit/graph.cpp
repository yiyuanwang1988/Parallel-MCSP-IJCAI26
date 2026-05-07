#include "graph.hpp"

constexpr int BITS_PER_UNSIGNED_INT (CHAR_BIT * sizeof(unsigned int));

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

rrsGraph::Graph::Graph(unsigned int n) {
    this->n = n;
    label = std::vector<unsigned int>(n, 0u);
    adjmat = {n, std::vector<unsigned int>(n, false)};
    degree = std::vector<unsigned int>(n, 0u);
    adjlist = {n, std::vector<unsigned int>(n, false)};
}

rrsGraph::Graph rrsGraph::induced_subGraph(struct rrsGraph::Graph& g, std::vector<int> vv) {
    Graph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<subg.n; j++)
            subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];

    for (int i=0; i<subg.n; i++)
        subg.label[i] = g.label[vv[i]];
    return subg;
}

void add_edge(rrsGraph::Graph& g, int v, int w, bool directed=false, unsigned int val=1) {
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
        //g.label[v] |= (1u << (BITS_PER_UNSIGNED_INT-1));
    }
}

static struct rrsGraph::Graph readDimacsGraph(char* filename, bool directed, bool vertex_labelled) {
    struct rrsGraph::Graph g(0);

    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL){
        fail("Cannot open file");
    }

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
                g = rrsGraph::Graph(nvertices);
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

static struct rrsGraph::Graph readLadGraph(char* filename, bool directed) {
    struct rrsGraph::Graph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int w;

    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = rrsGraph::Graph(nvertices);

    for (int i=0; i<nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)
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

static int read_word(FILE *fp) {
    unsigned char a[2];
    if (fread(a, 1, 2, fp) != 2)
        fail("Error reading file.\n");
    return (int)a[0] | (((int)a[1]) << 8);
}

static struct rrsGraph::Graph readBinaryGraph(char* filename, bool directed, bool edge_labelled,
        bool vertex_labelled)
{
    
    struct rrsGraph::Graph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "rb"))==NULL)
        fail("Cannot open file");

    int nvertices = read_word(f);
    g = rrsGraph::Graph(nvertices);

    // Labelling scheme: see
    // https://github.com/ciaranm/cp2016-max-common-connected-subGraph-paper/blob/master/code/solve_max_common_subGraph.cc
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

struct rrsGraph::Graph rrsGraph::readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled) {
    struct rrsGraph::Graph g(0);
    if (format=='D') g = readDimacsGraph(filename, directed, vertex_labelled);
    else if (format=='L') g = readLadGraph(filename, directed);
    else if (format=='B') g = readBinaryGraph(filename, directed, edge_labelled, vertex_labelled);
    else fail("Unknown Graph format\n");
    return g;
}

// void rrsGraph::set_adjlist(struct rrsGraph::Graph & g){
//     g.degree = new unsigned int [g.n];
//     g.adjlist = new unsigned int *[g.n];
//     unsigned int count = 0;
//     for(int i=0;i<g.n;++i){
//         count = 0;
//         for(int j=0;j<g.n;++j) if(g.adjmat[i][j]==1) count++;
//         g.adjlist[i] = new unsigned int[count]; count = 0;
        
//         for(int j=0;j<g.n;++j) 
//             if(g.adjmat[i][j]==1){
//                 g.adjlist[i][count]=j;
//                 ++count;
//             }
//         g.degree[i] = count; 
//     }
// }

void rrsGraph::set_adjlist(struct rrsGraph::Graph & g) {
    g.adjlist.clear();
    g.degree.clear();
    g.adjlist.resize(g.n);
    g.degree.resize(g.n);

    for (int i = 0; i < g.n; ++i) {
        g.adjlist[i].clear();
        for (int j = 0; j < g.n; ++j) {
            if (g.adjmat[i][j] == 1) {
                g.adjlist[i].push_back(j); 
            }
        }
        g.degree[i] = static_cast<ui>(g.adjlist[i].size());
    }
}

// void rrsGraph::GetEqClass(rrsGraph::Graph & g, ui *&EqClass){
//     ui Graph_size = g.n, label = 1, node = 0, node_neg = 0, node_nneg; EqClass = new ui[Graph_size];
//     bool equiv = true;
//     for(ui i=0;i<Graph_size;++i) EqClass[i]=0;
//     for(ui i=0;i<Graph_size;++i){
//         if(EqClass[i]!=0||g.degree[i]==0) continue;
//         node = g.adjlist[i][0];
//         for(ui j=0;j<g.degree[node];++j){
//             node_neg = g.adjlist[node][j];
//             if(g.degree[i]==g.degree[node_neg]&&node_neg>i){
//                 equiv=true;
//                 for(ui k=0;k<g.degree[node_neg];++k){
//                     node_nneg=g.adjlist[node_neg][k];
//                     if(node_nneg!=i&&g.adjmat[i][node_nneg]==0){
//                         equiv = false; break;
//                     }
//                 }
//                 if(equiv) EqClass[node_neg] = label;
//             }
//         }
//         EqClass[i]=label;
//         ++label;
//     }
// }


void rrsGraph::GetEqClass(rrsGraph::Graph & g, std::vector<int> &eq_class) {
    ui graph_size = g.n;
    eq_class.clear();
    eq_class.resize(graph_size, 0);
    if (graph_size == 0) return;

    ui label = 1;
    for (ui i = 0; i < graph_size; ++i) {
        if (eq_class[i] != 0 || g.degree[i] == 0) continue;
        if (g.adjlist[i].empty()) continue;
        ui node = g.adjlist[i][0];
        for (ui j = 0; j < g.degree[node]; ++j) {
            ui node_neg = g.adjlist[node][j];
            if (g.degree[i] == g.degree[node_neg] && node_neg > i) {
                bool equiv = true;
                for (ui k = 0; k < g.degree[node_neg]; ++k) {
                    ui node_nneg = g.adjlist[node_neg][k];
                    if (node_nneg != i && g.adjmat[i][node_nneg] == 0) {
                        equiv = false;
                        break;
                    }
                }
                if (equiv) eq_class[node_neg] = label;
            }
        }
        eq_class[i] = label;
        ++label;
    }
}
