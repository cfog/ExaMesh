#include <iostream>
#include <cstdio>

void 
inline writeAllTimeResults (FILE* file, int nP, double Tpartition, 
double TpartFaceMatching, double Textraction, double Trefinement, 
double TtriMatch, double TquadMatch, double Ttotal)
{
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    if (size == 0) 
	{
	 	fprintf(file, "%-5s %-12s %-18s %-12s %-12s %-12s %-12s %-12s\n",
          "nP", "Tpartition", "TpartFaceMatch", 
         "Textract","Trefine","TtriMatch","TquadMatch","Ttoal");
   	}
	
    fprintf(file, "%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f\n",
    nP,Tpartition,TpartFaceMatching,Textraction,Trefinement, 
    TtriMatch,TquadMatch,Ttotal);
}