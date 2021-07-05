#define LLEN 300
#define NAME 100
#define MAX_M 5000
#define SQR(x) (x)*(x)
typedef double vec_mt[3];
typedef struct
{
        int tag;
        vec_mt oxygen;
        vec_mt hydrogen1;
        vec_mt hydrogen2;
        int h_or_not;
} water_mt;
typedef struct
{
        int yes_or_no;
        int nodeid[1000];
        int node_number;
} cluster_mt;
double mindist(vec_mt point1, vec_mt point2, vec_mt boxlength)
{
        int i, j;
        vec_mt point3;
        for(i=0; i<3; i++)
                if(point1[i]-point2[i]>0.5*boxlength[i])
                        point3[i] = point1[i]-boxlength[i];
                else if(point1[i]-point2[i]<-0.5*boxlength[i])
                        point3[i] = point1[i]+boxlength[i];
                else
                        point3[i]=point1[i];
        return sqrt(SQR(point2[0]-point3[0])+SQR(point2[1]-point3[1])+SQR(point2[2]-point3[2]));
}
