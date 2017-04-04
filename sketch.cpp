#include <iostream>
#include <string>
#include <assert.h>
#include <vector>
#include <time.h>
#include <stdio.h>
#define DEBUGG 1
using namespace std;
const int INFINITY = -1; /*approximation for positive infinity*/

class Vertex
{
    public:
        int vertex;
        int distance;
        Vertex(int nVertex, int nDistance)
        {
            vertex=nVertex;
            distance=nDistance;
        }
    friend bool operator < (const Vertex & v1, const Vertex & v2)
    {
        if(v1.distance==INFINITY && v2.distance >=0)
        {
            return false;
        }
        else if (v1.distance==INFINITY && v2.distance==INFINITY)
        {
            return v1.vertex<v2.vertex;
        }
        else if (v1.distance>=0 && v2.distance==INFINITY)
        {
            return true;
        }
        else
        {
            if(v1.distance<v2.distance)
            {
                return true;
            }
            else if (v1.distance> v2.distance)
            {
                return false;
            }
            else
            {
                return v1.vertex<v2.vertex;
            }
        }

    }
    friend bool operator >(const Vertex & v1, const Vertex & v2)
    {
        if(v1.distance==INFINITY && v2.distance >=0)
        {
            return true;
        }
        else if (v1.distance==INFINITY && v2.distance==INFINITY)
        {
            return v1.vertex>v2.vertex;
        }
        else if (v1.distance>=0 && v2.distance==INFINITY)
        {
            return false;
        }
        else
        {
            if(v1.distance>v2.distance)
            {
                return true;
            }
            else if (v1.distance< v2.distance)
            {
                return false;
            }
            else
            {
                return v1.vertex>v2.vertex;
            }
        }

    }
};

/* this is intended for implememtation of priority queue supproting a log(n) time of operations including
extract min, decrease key for Dijkstra algorithm*/
/* this priority queue is build specifically for implementation of Dijkstra algorithm*/
/*for even fatser implementation, a fibonacci heap is an option*/
class BinaryHeap
{
    vector<Vertex> heap;
    int* vertexIndex; /* this array is aimed to maintained the indices of vertice in the heap, to support key decreasing*/
    private:
        int parent(int i) /* return the index of parent of i */
        {
            return i/2;
        }
        int left(int i) /*index of left child of i*/
        {
            return 2*i;
        }
        int right(int i) /*index of right child of i*/
        {
            return 2*i+1;
        }
        void swap(int i, int j) /*swap 2 elements in the heap, meanwhile modify the vertexIndex*/
        {
            /*modify vertexIndex first*/
            vertexIndex[heap[i].vertex]=j;
            vertexIndex[heap[j].vertex]=i;
            /*swap the elemnts*/
            Vertex temp=heap[i];
            heap[i]=heap[j];
            heap[j]=temp;
        }
        void minHeapify (int i)/*i is the node to be heapified, subtree assumed to obey the property*/
        {
            int l=left(i);
            int r=right(i);
            int smallest=i;
            if(l<heap.size() && heap[l]<heap[smallest])
            {
                smallest=l;
            }
            if(r<heap.size() && heap[r]<heap[smallest])
            {
                smallest=r;
            }
            if(smallest!=i)
            {
                swap(i, smallest);
                minHeapify(smallest);
            }
        }
        void decreaseWithIndex(int i, int nDistance)
        {
            Vertex temp(heap[i].vertex, nDistance);/*since the distance comparison involves tedious check of INFINITY, I will use temp vertex here*/
            if(temp < heap[i])
            {
                heap[i].distance=nDistance;
                while(i>0  &&  heap[parent(i)] > heap[i]) /*maintain the heap property*/
                {
                    swap(i, parent(i));
                    i=parent(i);
                }
            }
        }

    public:
        BinaryHeap(int n) /*n is the number of vertices*/
        {
            vertexIndex = new int[n];
            for(int i =0; i<n; ++i)
            {
                vertexIndex[i]=-1; /*all indices are set to be -1*/
            }
        }
        void deleteArray() /* to free memory allocated for the array vertexIndex */
        {
            delete[] vertexIndex;
        }
        bool empty()
        {
            return heap.size()==0;
        }
        bool contains (int vertex) /*test is the target vertex is contained in the heap*/
        {
            return vertexIndex[vertex]!=-1;/* when it is -1, then it is not contained in the heap*/
        }
        Vertex min() /*return the minimum, but not deleting the minimum*/
        {
            return heap[0];
        }
        Vertex extractMin()/*return and remove the minimum*/
        {
            Vertex min=heap[0];
            swap(0, heap.size()-1); /*swap the first and the last element in the heap*/
            vertexIndex[heap[heap.size()-1].vertex]=-1;/*the minimum is removed form the heap, then its index is -1 to avoid potential conflicts*/
            heap.pop_back(); /*minimum is now poped out of the heap */
            minHeapify(0); /*to ensure the heap property*/
            return min;
        }
        void decreaseWithVertex(Vertex targetVertex)
        {
            decreaseWithIndex(vertexIndex[targetVertex.vertex], targetVertex.distance);
        }
        void insert (Vertex vertexInsert) /*decreaseWithVertex function is used to ease the implementation*/
        {
            Vertex infinityVertex(vertexInsert.vertex, INFINITY);
            heap.push_back(infinityVertex);
            vertexIndex[vertexInsert.vertex]=heap.size()-1;
            decreaseWithVertex(vertexInsert);
        }
        void deleteVertexIndex()
        {
            delete[] vertexIndex;
        }
        int distanceOfVertex(int vertexID)
        {
            return heap[vertexIndex[vertexID]].distance;
        }
};



class GraphVertex {
    public:
        int vertex;
        int weight; /*this weight represent its distance to its neighbor, where neighbor is the indice in the array in the adjcent list structure*/
        GraphVertex* next;
        GraphVertex(int v, int w)
        {
            vertex=v;
            weight=w;
        }

};

class Graph
{
    private:
        GraphVertex** graph;
        int size;
        int* distance;
    public:
        Graph(int n)
        {
            size=n;
            graph=new GraphVertex*[size](); /*to initialize all the elements in the array as 0/NULL*/
            distance =new int[size];
            for(int i=0; i<size; ++i)
            {
                distance[i]=INFINITY;
            }
        }
        void addEdge(int v1, int v2, int w)
        {
            GraphVertex* temp1=new GraphVertex(v2, w);
            temp1->next=graph[v1];
            graph[v1]=temp1;
            GraphVertex* temp2=new GraphVertex(v1,w);
            temp2->next=graph[v2];
            graph[v2]=temp2;
        }
        void deleteGraph()
        {
            for(int i=0; i<size; ++i)
            {
                GraphVertex* current=graph[i];
                while(current!=NULL)
                {
                    GraphVertex* temp=current->next;
                    delete current;
                    current=temp;
                }
            }
            delete[] graph;
            delete[] distance;
        }
        void shortestPath(int sourceID)
        {
            /*given a vertex, finding the shortest pathes between it and all the other vertices, infinity if not connected*/
            /*path is not reserved, however, adding an array storing the previous node in the path would do*/
            int* visited = new int[size]();
            BinaryHeap priorityQ(size);
            distance[sourceID]=0;
            Vertex source(sourceID, 0);
            priorityQ.insert(source);
            while(priorityQ.empty()==false)
            {
                source=priorityQ.extractMin();
                sourceID=source.vertex;
                visited[sourceID]=1;
                GraphVertex* neighbor=graph[sourceID];
                while(neighbor!=NULL)
                {
                    int vertexID=neighbor->vertex;
                    if(visited[vertexID]==0)
                    {
                        /*for all vertices that is unvisited, if it is not in the priorityQ, its distance must be INFINITY, else it must not be INFINITY*/
                        int nDistance=distance[sourceID]+neighbor->weight;
                        Vertex nVertex(vertexID, nDistance);
                        if(priorityQ.contains(vertexID)==false)
                        {
                            distance[vertexID]=nDistance;
                            priorityQ.insert(nVertex);
                        }
                        else
                        {
                            if(distance[vertexID]>nDistance)
                            {
                                distance[vertexID]=nDistance;
                                priorityQ.decreaseWithVertex(nVertex);
                            }
                        }
                    }
                    neighbor=neighbor->next;
                }/*end of inserting/decreasing keys for the neighbor of source*/
            }/*end of looping for the priorityQ*/
            priorityQ.deleteVertexIndex(); /*if time is intense, remove this delete statement*/
        }/*end of the dijkstra algorithm*/
    
    
        int shortestPathBetween2points(int sourceID, int destinationID)
        {
            int* visited = new int[size]();
            BinaryHeap priorityQ(size);
            Vertex source(sourceID, 0);
            priorityQ.insert(source);
            while(priorityQ.empty()==false)
            {
                source=priorityQ.extractMin();
                sourceID=source.vertex;
                if(sourceID==destinationID)
                {
                    delete[] visited;/*if time is limited, remove this*/
                    return source.distance;
                }
                visited[sourceID]=1;
                GraphVertex* neighbor=graph[sourceID];
                while(neighbor!=NULL)
                {
                    int vertexID=neighbor->vertex;
                    if(visited[vertexID]==0)
                    {
                        /*for all vertices that is unvisited, if it is not in the priorityQ, its distance must be INFINITY, else it must not be INFINITY*/
                        /*the new distance can be calculated from the weight between it and the source vertex plus the source distance*/
                        int nDistance=source.distance+neighbor->weight;
                        Vertex nVertex(vertexID, nDistance);
                        if(priorityQ.contains(vertexID)==false)
                        {
                            priorityQ.insert(nVertex);
                        }
                        else
                        {
                            if(priorityQ.distanceOfVertex(vertexID)>nDistance)
                            {
                                priorityQ.decreaseWithVertex(nVertex);
                            }
                        }
                    }
                    neighbor=neighbor->next;
                }/*end of inserting/decreasing keys for the neighbor of source*/
            }/*end of looping for the priorityQ*/
            priorityQ.deleteVertexIndex();/*if time is limited, remove this*/
            delete[] visited;/*if time is limited, remove this*/
            return INFINITY;
        }/*this is a littel amendment of dijkstra algorithm*/
    
        int distanceOfVertex(int u)
        {
            return distance[u];
        }

};

void NearestDriver(){
    int n, m;

    cin >> n >> m;

    Graph g(n); // implement Graph class by yourself

    for(int i = 0; i < m; i++){
        int a, b, w;
        cin >> a >> b >> w;

        g.addEdge(a, b, w);
    }

    int u;
    cin >> u;
    // implement your own shortest path
    g.shortestPath(u);

    int bestv = -1;
    int l;
    cin >> l;
    for(int i = 0; i < l; i++){
        // scan over every car to get the final answer
        int v;
        cin>>v;
        int vdistance=g.distanceOfVertex(v);
        if(vdistance!=INFINITY)
        {
            if(bestv==-1)
            {
                bestv=v;
            }
            else
            {
                int bdistance=g.distanceOfVertex(bestv);
                if(bdistance==INFINITY)
                {
                    bestv=v;
                }
                else if(bdistance==vdistance)
                {
                    if(bestv>v)
                    {
                        bestv=v;
                    }
                }
                else if(bdistance>vdistance)
                {
                    bestv=v;
                }
            }
        }
    }

    if(bestv == -1)
        cout << "NO" << endl;
    else
        cout << bestv << endl;
    g.deleteGraph();/*if time is intense, remove this delete statement*/
    return;
}


void QueryPrice(){
    //your code starts here
    int n, m;
    cin >> n >> m;
    
    Graph g(n); // implement Graph class by yourself
    
    for(int i = 0; i < m; i++){
        int a, b, w;
        cin >> a >> b >> w;
        
        g.addEdge(a, b, w);
    }
    int k;
    cin>>k;
    for(int i=0; i<k; ++i)
    {
        int s, t;
        cin>>s>>t;
        int distance=g.shortestPathBetween2points(s, t);
        if(distance!=INFINITY)
        {
            cout<<distance<<endl;
        }
        else
        {
            cout<<"NO"<<endl;
        }
    }
    g.deleteGraph();/*if time is limited, remove this delete statement*/
    return;
}

void Diameter(){
    //your code starts here
}

void DiameterApproximation(){
    //your code starts here
}

int main(){
    #if DEBUGG
        clock_t tStart = clock();
    #endif
    string section;
    cin >> section;

    if(section == "NearestDriver")
        NearestDriver();
    else if(section == "QueryPrice")
        QueryPrice();
    else if(section == "Diameter")
        Diameter();
    else if(section == "DiameterApproximation")
        DiameterApproximation();
    else{
        cout << "wrong input file!" << endl;
        assert(0);
    }
    #if DEBUGG
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    #endif
    return 0;
}

