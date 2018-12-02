    #include <stdio.h>
    #include <dirent.h>
    #include <string>
    #include <iostream>
    #include <vector>
    #include <fstream>
    #include <cstring>
    #include <map>
    #include <stdlib.h>
    #include <algorithm>
    #include "Shingles.hpp"
	#include <set>
	#include <queue>
    #include <iomanip>
 
    


    using namespace std;
    typedef pair< int, int > pii;
    typedef pair< string , pii > p_str_pii;
    typedef vector <string> vs;
    typedef vector< vs > Sequences;
    using pfi = pair<float,int>;
    using pif = pair <int,float>;
    using vi = vector<int>;
    using msi = map<string,int>;
    using mis = map<int,string>;
    using mii = map<int,int>;
    using vf = vector < float >;


    
    vs obtener_nombres_archivos(int c, char *v[]);
    vs split(string line, string delims);
    vs archivos_validos(vs aux);
    mis correlacion(vs nombres_archivos);


    int main (int c, char *v[]) {

        vs archivos_en_directorio=obtener_nombres_archivos(c,v),complex,nombres_archivos=archivos_validos(archivos_en_directorio);
        vector<vs> complejos;
        map<string,vector<pfi>> recall_doc,precision_doc;
        vector< p_str_pii > lista;
        vf recall,precision;
        vector < vf > vrecall, vprecision;
        vector<vector <p_str_pii>> completo;
        map<vs,vi> mapita;   
        msi docs,map_complejos;
        mis correlativo=correlacion(nombres_archivos);
        mii documents,doc_no_asociados;

        bool encontro_complex_en_doc =false;
        int contador_linea=1,numdoc=0;
        // Se limpia los nombres de archivos filtrandolos por los que contengan substring ".txt"


        //Se itera sobre cada documento (Complejo Predicho)      

        for(string i : nombres_archivos){
            contador_linea=1;
       		encontro_complex_en_doc=false;
            string aux;
            int list_size;

            
            ifstream fe(i); // abre el archivo
            while(getline (fe,aux)) { // obtiene linea desde el texto
                vs vstr=split(aux,"\t,");                
                if(vstr.size()!=0 && vstr[0].compare("List")==0){ // Si llega a la septima linea guarda <nombre_doc,list_size> y <ID_doc,List_size>
                       list_size=stoi(vstr[1]);
                       docs.insert({i,list_size});
                       documents.insert({numdoc,list_size});
                }
                   
                if(vstr.size()!=0 && vstr[1].find("complex")!=string::npos){
 		
                        int hits=stoi(vstr[4]),real=stoi(vstr[5]);
                        string complejo = vstr[1];
                        float rec=(float)hits/(float)real,preci= (float)hits/(float)list_size;
                        pii nuevo = {hits,real}; //guarda List Annotation y List BG annotation
                        recall.push_back(rec);
                        precision.push_back(preci);
 						lista.emplace_back(complejo,nuevo);
 						encontro_complex_en_doc=true;
 						complex.push_back(complejo); // se almacena un vector conteniendo a todos los complejos encontrados
 						map_complejos.insert({complejo,real}); // se mapea el nombre del complejo real a la cantidad total de proteinas
 						if(recall_doc.find(complejo)!=recall_doc.end()){
     						recall_doc[complejo].push_back(pfi(rec,numdoc)); 
                            precision_doc[complejo].push_back(pfi(preci,numdoc));
 						}else{
     						vector<pfi> plop,plap;
     						plop.push_back(pfi(rec,numdoc));
                            plap.push_back(pfi(preci,numdoc));
     						recall_doc[complejo]=plop;
                            precision_doc[complejo]=plap;
 						}
 					}
                
                contador_linea++;
            }
            if (!encontro_complex_en_doc) // en las siguientes lineas la informacion recuperada antes es insertada a sus vectores correspondientes
            {
            	pair<int,int> n=make_pair(-1,-1);
            	lista.emplace_back("none",n);
            	completo.push_back(lista);
            	vs aux ;
            	aux.push_back("none");
            	//complejos.push_back(aux);
                //agregar los complejos que no tiene ninguna coincidencia con algun complejo real a un mapa
                doc_no_asociados[numdoc]=list_size;
                if(mapita.find(aux)!=mapita.end()){
                   
                    mapita[aux].push_back(numdoc);

                 }else{
                    vector <int> vectorcito;
                    vectorcito.push_back(numdoc);
                    mapita.insert({aux,vectorcito});
                }
            }else {
             completo.push_back(lista);          
             complejos.push_back(complex);
             vrecall.push_back(recall);
             vprecision.push_back(precision);
             if(mapita.find(complex)!=mapita.end()){
                vi vec=mapita[complex];
                vec.push_back(numdoc);
                mapita.insert({complex,vec});
             }else{
                vector <int> vectorcito;
                vectorcito.push_back(numdoc);
                mapita.insert({complex,vectorcito});
            }
             
         	}
            lista.clear();
            complex.clear();
            fe.close();
            numdoc++;
        }
      
        int n=0;
        ofstream fx("complejos_real_mapeado_y_recall");
        for(const auto& aux : recall_doc){

            fx<<aux.first<<endl;
            for(const auto& aux2: aux.second) fx<<"ID_doc: "<<aux2.second<<'\t'<<"Recall: "<<aux2.first<<endl;
            fx<<endl;
        }
        fx.close();
        Sequences sequences = complejos;

        // Build clusters
        map<Shingles::Signature, Sequences> clusters;
        Shingles shingle(3);

        for (const auto& seq : sequences) {
           
        Shingles::Signature signature = shingle.compute(seq);

        clusters[signature].push_back(seq);
        
         }

        // Guardamos los Cluster en un archivo de salida
         ofstream fs("Clusters");    
        for (const auto& cluster : clusters) {
            for (const auto& seq : cluster.second) {
                fs << cluster.first << ": ";
                for (const auto& val : seq)
                   fs << val << " ";
                fs << '\n';
            }
        fs << '\n';
        n++;
        }
        fs.close();
        map<vs,pif> mejor_por_seq;

        ofstream fp("Mejor por secuencia");
         vector<priority_queue<pfi>> cola1;
        // vector<priority_queue<pfi>> cola2;
          for (const auto& cluster:clusters)
            {   priority_queue<int > aux;
                priority_queue <pfi> colita1;
                priority_queue <pfi> colita2;

                for(const auto& seq: cluster.second){
                    fp<< cluster.first <<":"<<'\t';
                    vi help= mapita[seq];
                    pair < int , int > max=make_pair(0,0);
                    float maximum=-1;
                    int docum;
                    for(int i:help){

                        for(const auto& j:completo[i]){
                             if(((float)j.second.first/(float)j.second.second)>(float)maximum){

                                    max=j.second; 
                                    docum=i;
                                    maximum=((float)j.second.first/(float)j.second.second);
                                }
                        }
                     /*   fp<<"complejo :"<<j[0].first<<"hits:"<<j[0].second.first<<"    total real:"<<j[0].second.second<<endl;
                        float precision=(float)j[0].second.first/(float)documents[i];
                        float recall=(float)j[0].second.first/(float)j[0].second.second;

                        fp<<"Precision :"<<setw(9)<<precision<<"   Recall:"<<setw(9)<<recall<<endl;
                        colita1.push({precision,i});
                        colita2.push({recall,i});*/
                    }
                mejor_por_seq[seq]=make_pair(maximum,docum);  
                colita1.push(make_pair(maximum,docum));  
                fp<<"mejor_doc_por_el_momento_dentro_de_cluster: "<<colita1.top().second<<"    recall: "<<colita1.top().first<<endl;    
            }
            cola1.push_back(colita1);
            //cola2.push_back(colita2);
    }
    int counter=0;
    for (const auto& cluster : clusters) {
                fp<<left<<setw(20)<<cluster.first<<": Documento:"<<cola1[counter].top().second<<"   Recall:"<<cola1[counter].top().first<<endl;
                counter++;
    }

    fp.close();
    //crearemos un vector de entero rellenado con los numeros de documentos
    vector <int> vi;
    int jo=0;
    for(const auto& j: nombres_archivos){
        vi.push_back(jo);
        jo++;
    }

    map<string,pif> final;

    map<string,bool> checkeado;

 ofstream fn("Mejor por complejo real");
  
  for(const auto& cluster:clusters){
        for(const auto& seq:cluster.second){
            for(const auto& val : seq){
                if(checkeado.find(val)==checkeado.end()){

                vector < pfi> var = recall_doc[val];
                sort(var.begin(), var.end(), [](const pfi &left, const pfi &right) {
                return left.first < right.first;});
                fn <<left <<setw(65)<<val<<"  Recall:"<<setw(15)<<var.back().first <<"    Documento:"<<var.back().second<<endl;
                final[val]=make_pair(var.back().first,var.back().second);
                //recall_doc[val].erase(remove(recall_doc[val].begin(),recall_doc[val].end(),var.back()), recall_doc[val].end());
                checkeado[val]=true;
                fn<<endl;
                }
            }
        }            
    }
      fn.close();
  /*  for (const auto& fin:final)
    {
       // fn<<"Complejo: "<<left<<setw(65)<<fin.first<<" "<<"documento : "<<fin.second.second << "     Recall: "<<fin.second.first<<endl;
    }*/
    ofstream p1("Complejos no asociados");
    for(const auto& i:doc_no_asociados){
        p1<<"Complejo: "<<i.first<<"        archivo: "<<left<<setw(60)<<correlativo[i.first]<<" List size: "<<i.second<<endl;

    }
    p1<<"Cantidad: "<<doc_no_asociados.size()<<endl;
    p1.close();
    return 0;
}

    vs obtener_nombres_archivos(int c, char *v[]){

        vs list_archivos;
        struct dirent *pDirent;
        DIR *pDir;
        string str = "";

        if (c < 2) {
            printf ("Usage: testprog <dirname>\n");
            return list_archivos;
        }
        pDir = opendir (v[1]);
        if (pDir == NULL) {
            printf ("Cannot open directory '%s'\n", v[1]);
            return list_archivos;
        }

        while ((pDirent = readdir(pDir)) != NULL) {
           list_archivos.push_back(pDirent->d_name);

        }
        closedir (pDir);

        return list_archivos;
    }

    vs split(string line, string delims){
       string::size_type bi, ei;
       vs words;
       bi = line.find_first_not_of(delims);
       while(bi != string::npos) {
            ei = line.find_first_of(delims, bi);
            if(ei == string::npos)
                    ei = line.length();
                    words.push_back(line.substr(bi, ei-bi));
                    bi = line.find_first_not_of(delims, ei);
            }
      return words;
    }
    vs archivos_validos(vs aux){
        vs nombres_archivos;
        for ( string i : aux){
            if(i.find(".txt")!=string::npos){
                nombres_archivos.push_back(i);
            }
        }
    return nombres_archivos;
    }
    mis correlacion(vs nombres_archivos){
        mis salida;
        int max=nombres_archivos.size();
        for(int i=0;i<max;i++){
            salida[i]=nombres_archivos[i];
        }

        return salida;
    }
