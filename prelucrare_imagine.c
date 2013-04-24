/*
 * Voiculescu Rafaela
 * 331CA
 *
 * APD - Tema 3 - Prelucrare de imagini
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

#define DIM 100 //dimensiunea unei linii citite dintr-un fisier

#define DATE 1  //tag pentru fasia primita (liniile din matricea img initiala) de procesul p
#define SOL 2   //tag pentru trimiterea liniilor dupa prelucrarea realizatata e fiecare proces
#define DATE_FILTRU 3   //tag pentru datele schimbate intre ele de procese pentru a realiza filtrul asupra imaginii
#define MAX_POSIBIL 32000
#define MIN_POSIBIL -32000

int nrProc;         //numar de procese pornite
FILE *fin,*fout;    //fisiere de intrare, respectiv iesire
int **img;          //matricea pentru imagine
int col,lin,valMax; //nr de coloane, respectiv linii al matricei si valoarea maxima ce apare in ea
int filtru[][9]={{1,1,1,1,1,1,1,1,1},
                 {1,2,1,2,4,2,1,2,1},
                 {0,-2,0,-2,11,-2,0,-2,0},
                 {-1,-1,-1,-1,9,-1,-1,-1,-1},
                 {-1,0,-1,0,4,0,-1,0,-1}};
int impartit[]={9,16,3,1,1};
int adaugare[]={0,0,0,0,127};

//functie pentru citirea din fisierul de intrare a datelor
void citireDateIntrare()
{
    int ok=0,i,j;
    long pozCrt;
    long offset;
    //int col,lin,valMax;
    char linie[DIM];

    //se face citire pana se ajunge la o linie ce incepe cu P sau se termina fisierul
    while(fgets(linie,DIM+1,fin)!=NULL)
    {
        //debug
        //printf("!%s\n",linie);

        if(linie[0]=='P')   //s-a gasit linia care incepe cu P
        {
            ok=1;
            break;
        }
    }

    if(ok==0)   //s-a ajus la finalul fisierului si nu s-a gasit o linie cu P
    {
        printf(">>>>>>>>>>Eroare: In fisier nu apare linia care incpe cu litera 'P'!<<<<<<<<<<");
        MPI_Finalize();
        exit(-1);
    }
    
    //verificare daca poza este P2
    //printf("DIM linie: %d \n",strlen(linie));
    if(linie[1]!='2')
    {
        //debug
        //puts("Afisare linie:");
        //for(i=0;i<3;i++)
        //   printf("%d ",(int)linie[i]);
        //printf("!%c\n",linie[strlen(linie)-1]);

        printf(">>>>>>>>>>Eroare: Pe prima linie nu se afla P2!<<<<<<<<<<");
        MPI_Abort(MPI_COMM_WORLD,-1);
        exit (-1);
    }

    //debug
    fprintf(fout,"%s",linie);


    //se sar liniile care nu incep cu # sau cifra (liniile goale)
    while(fgets(linie,DIM+1,fin)!=NULL && linie[0]!='#' && (linie[0]>'9' || linie[0]<'0'))
    {
        //debug
        //printf("!%s\n",linie);
    }

    //verificam daca linia incepe cu cifra -> ne intoarcem la inceputul ei si abia apoi citim
    if(linie[0]<'9' && linie[0]>'0')
    {
        pozCrt=ftell(fin);             //aflam pozitia curenta a cursorului in fisier
        //debug
            //printf("Dim linie cu cifre %d\n",strlen(linie));
        //offset=pozCrt-strlen(linie)-1; //pozitia inceputului liniei curente in fisier
        offset=pozCrt-strlen(linie);        //pozitia inceputului liniei curente in fisier
        fseek(fin,offset,SEEK_SET);         //deplasare cursor la inceputul liniei curent in fisier
    }

    //se verifica daca linia incepe cu #
    //if(linie[0]=='#')
    //linia incepe cu # -> citim mai departe
    {

        //debug
        //if(linie[0]=='#')
            //printf("Dim linie cu # %d\n",strlen(linie));
        //aflam pozitia din fisier
        //pozCrt=ftell(fin);

        //debug
        //if(linie[0]=='#')
           // printf("poz cursor in fis cu #: %d\n",(int) pozCrt);
        if(linie[0]=='#')
            fprintf(fout,"%s",linie);
        //printf("%s",linie);

        //pozitionam cursorul inainte de numar
        //fseek(fin,pozCrt-1,SEEK_SET);

        //citim nr coloane al matricei
        fscanf(fin,"%d",&col);

        //citim nr linii al matricei
        fscanf(fin,"%d",&lin);

        //citim valoarea maxima ce va aparea in matrice
        fscanf(fin,"%d",&valMax);

        //debug
        //printf("Nr col, nr lin, valMax: %d %d %d \n",col,lin,valMax);
        
        fprintf(fout,"%d %d\n%d",col,lin,valMax);

        //alocam spatiu pentru matrice
        img=(int **) calloc (lin+2,sizeof(int*));
        for(i=0;i<lin+2;i++)
            img[i]=(int *) calloc (col+2,sizeof(int));

        //citim matricea
        for(i=1;i<=lin;i++)
            for(j=1;j<=col;j++)
            {
                fscanf(fin,"%d",&img[i][j]);
                //printf("%d ",img[i][j]);
            }
        //debug
        fprintf(fout,"%s","\n");
        //for(i=1;i<=lin;i++)
        //{
        //    for(j=1;j<=col;j++)
        //        fprintf(fout,"%d\n",img[i][j]);
        //}

        //return;
    }   
}

void contrastImagine(int rank,int nrLin,int nrCol,int **matrice, int a,int b)
{
    //debug
    //printf("------------START prelucrare imagine pentru procesul %d----------------\n",rank);

    //debug
    //printf("a,b: %d %d\n",a,b);

    int minLoc=MAX_POSIBIL; //va retine minimul din fasia prelucrata de procesul rank
    int maxLoc=MIN_POSIBIL; //va retine maximul din fasia prelucrata de procesul rank
    int min;    //va retine minimul global (din toata poza)
    int max;    //va retine maximul global (din toata poza)

    int i,j;

    for(i=1;i<=nrLin-2;i++)
        for(j=1;j<=nrCol-2;j++)
        {
            if(matrice[i][j]>maxLoc)
            {
                maxLoc=matrice[i][j];
            }
            if(matrice[i][j]<minLoc)
            {
                minLoc=matrice[i][j];
            }
        }
    //debug
    //printf("***Minim/Maxim pt proces %d : %d/%d\n",rank,minLoc,maxLoc);

    //punem bariera sa fim siguri ca fiecare proces si-a calculat maximul respectiv minimul global
    MPI_Barrier(MPI_COMM_WORLD);

    //aflam minim global
    MPI_Allreduce(&minLoc,&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    //aflam maxim global
    MPI_Allreduce(&maxLoc,&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    //debug
    //printf(">>>Minim/Maxim GLOBALE calculate pt proces %d : %d/%d\n",rank,min,max);

    //parcurgem matricea corespunzatoare procesului si facem modificarile necesare
    for(i=1;i<=nrLin-2;i++)
        for(j=1;j<=nrCol-2;j++)
            matrice[i][j]=(b-a)*(matrice[i][j]-min)/(max-min)+a;

    //debug
   // printf("------------END prelucrare imagine pentru procesul %d----------------\n",rank);
}

void filtruImagine(int rank,int nrLin,int nrCol,int **matrice,char* fil)
{
    MPI_Status stat;
    int i,j;
    int sum=0;
    int tip;    //tipul de filtru aplicat
    int **matr_aux;
    
    //comunicarea intre procese pentru a isi completa datele curente

    //toate procesele in afara de procesul 0 primesc de la procesul anterior date
    if(rank!=0)
    {
        MPI_Recv(matrice[0],nrCol,MPI_INT,rank-1,DATE_FILTRU,MPI_COMM_WORLD,&stat);
    }
    //toate procesele in afara de ultimul trimit date procesului ce le succede
    if(rank!=nrProc-1)
    {
        MPI_Send(matrice[nrLin-2],nrCol,MPI_INT,rank+1,DATE_FILTRU,MPI_COMM_WORLD);
    }
    //toate procesele in afara de ultimul primesc date de la cel de dupa ele
    if(rank!=nrProc-1)
    {
        MPI_Recv(matrice[nrLin-1],nrCol,MPI_INT,rank+1,DATE_FILTRU,MPI_COMM_WORLD,&stat);
    }
    //toate procesele in afara de primul trimit date celui de inaintea lor
    if(rank!=0)
    {
        MPI_Send(matrice[1],nrCol,MPI_INT,rank-1,DATE_FILTRU,MPI_COMM_WORLD);
    }

    //debug -afisare matrice completa
    /*printf("---> matricea completa pentru procesul %d <---\n",rank);
    for(i=0;i<nrLin;i++)
    {
        for(j=0;j<nrCol;j++)
            printf("%d ",matrice[i][j]);
        printf("\n");
    }
    printf("<--- matricea completa pentru procesul %d --->\n",rank);*/
    //alocarea memoriei pentru matricea auxiliara
    matr_aux=(int **) calloc(nrLin,sizeof(int*));
    for(i=0;i<nrLin;i++)
        matr_aux[i]=(int *) calloc(nrCol,sizeof(int));

    //verificam ce tip de filtru trebuie aplicat
    if(strcasecmp(fil,"smooth")==0)
        tip=0;
    else
        if(strcasecmp(fil,"blur")==0)
            tip=1;
        else
            if(strcasecmp(fil,"sharpen")==0)
                tip=2;
            else
                if(strcasecmp(fil,"mean_removal")==0)
                    tip=3;
                else
                    if(strcasecmp(fil,"emboss")==0)
                        tip=4;
                    else
                    {
                        printf(">>>>>>>>>> Filtrul specificat nu este corect <<<<<<<<<<");
                        MPI_Abort(MPI_COMM_WORLD,-1);
                        exit(-1);
                    }

    //aplicam filtrul pe matrice
    for(i=1;i<=nrLin-2;i++)
        for(j=1;j<=nrCol-2;j++)
        {
            sum=0;
            sum=sum+matrice[i-1][j-1]*filtru[tip][0];
            sum=sum+matrice[i-1][j]*filtru[tip][1];
            sum=sum+matrice[i-1][j+1]*filtru[tip][2];
            sum=sum+matrice[i][j-1]*filtru[tip][3];
            sum=sum+matrice[i][j]*filtru[tip][4];
            sum=sum+matrice[i][j+1]*filtru[tip][5];
            sum=sum+matrice[i+1][j-1]*filtru[tip][6];
            sum=sum+matrice[i+1][j]*filtru[tip][7];
            sum=sum+matrice[i+1][j+1]*filtru[tip][8];

            sum/=impartit[tip];

            sum+=adaugare[tip];

            //retinem noua valoare in matricea auxiliara
            if(sum<0)
                matr_aux[i][j]=0;
            else
                if(sum>valMax)
                    matr_aux[i][j]=valMax;
                else
                    matr_aux[i][j]=sum;

        }

    //facem update pentru informatiile din matrice dupa ce au fost aplicate modificarile
    for(i=0;i<nrLin;i++)
        for(j=0;j<nrCol;j++)
            matrice[i][j]=matr_aux[i][j];

    //dezalocarea memoriei pentru matricea auxiliara
    for(i=0;i<nrLin;i++)
        free(matr_aux[i]);
    free(matr_aux);
}

void entropieImagine(int rank,int nrLin,int nrCol,int **matrice, int a,int b,int c)
{
    MPI_Status stat;
    int i,j;
    int sum=0;
    int **matr_aux;

    //comunicarea intre procese pentru a isi completa datele curente

    //toate procesele in afara de procesul 0 primesc de la procesul anterior date
    if(rank!=0)
    {
        MPI_Recv(matrice[0],nrCol,MPI_INT,rank-1,DATE_FILTRU,MPI_COMM_WORLD,&stat);
    }
    //toate procesele in afara de ultimul trimit date procesului ce le succede
    if(rank!=nrProc-1)
    {
        MPI_Send(matrice[nrLin-2],nrCol,MPI_INT,rank+1,DATE_FILTRU,MPI_COMM_WORLD);
    }
    //toate procesele in afara de ultimul primesc date de la cel de dupa ele
    if(rank!=nrProc-1)
    {
        MPI_Recv(matrice[nrLin-1],nrCol,MPI_INT,rank+1,DATE_FILTRU,MPI_COMM_WORLD,&stat);
    }
    //toate procesele in afara de primul trimit date celui de inaintea lor
    if(rank!=0)
    {
        MPI_Send(matrice[1],nrCol,MPI_INT,rank-1,DATE_FILTRU,MPI_COMM_WORLD);
    }

    //debug -afisare matrice completa
    //printf("---> matricea completa pentru procesul %d <---\n",rank);
    /*for(i=0;i<nrLin;i++)
    {
        for(j=0;j<nrCol;j++)
            printf("%d ",matrice[i][j]);
        printf("\n");
    }*/
    //printf("<--- matricea completa pentru procesul %d --->\n",rank);
    //alocarea memoriei pentru matricea auxiliara
    matr_aux=(int **) calloc(nrLin,sizeof(int*));
    for(i=0;i<nrLin;i++)
        matr_aux[i]=(int *) calloc(nrCol,sizeof(int));

    //calcul entropie
    for(i=1;i<=nrLin-2;i++)
        for(j=1;j<=nrCol-2;j++)
            matr_aux[i][j]=matrice[i][j]-ceil(a*matrice[i-1][j] + b*matrice[i-1][j-1] + c*matrice[i][j-1]);

    //facem update pentru informatiile din matrice dupa ce au fost aplicate modificarile
    for(i=0;i<nrLin;i++)
        for(j=0;j<nrCol;j++)
            matrice[i][j]=matr_aux[i][j];

    //dezalocarea memoriei pentru matricea auxiliara
    for(i=0;i<nrLin;i++)
        free(matr_aux[i]);
    free(matr_aux);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);

    int rank;
    int i,j,k;
    int nrLin,nrCol;    //numarul de linii respectiv coloane al matricei fiecarui proces (fasie din cea initiala)
    int nrLinProc;

    MPI_Status status;

    int **matrice;

    //aflam numarul de procese pornite
    MPI_Comm_size(MPI_COMM_WORLD,&nrProc);
    //debug- afisare numar procese pornite
    //printf(">>>>>>>>>>Numar de procese pornite: %d<<<<<<<<<<\n",nrProc);

    //aflam rank-ul procesului curent
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank==0) //daca e proces master
    {
        //deschidem fisier e intrare in functie de tipul de ajusater pe care-l facem imaginii (ordinea parametriclor e alta)
        if(atoi(argv[1])==0)    //constrast
        {
            fin=fopen(argv[2],"r");
            fout=fopen(argv[5],"w");
        }
        else
            if(atoi(argv[1])==1)    //filtru
            {
                fin=fopen(argv[2],"r");
                fout=fopen(argv[4],"w");
            }
            else
                if(atoi(argv[1])==2)    //entropie
                {
                    fin=fopen(argv[2],"r");
                    fout=fopen(argv[6],"w");
                }
                else
                {
                    printf(">>>>>>>>>>Parametrul dat nu este corect! 0-contrast; 1-filtru; 2-entropie <<<<<<<<<<");
                    MPI_Abort(MPI_COMM_WORLD,-1);
                    exit(-1);                    
                }

        //verificare ca s-a deschis fisierul de intrare
        if(fin==NULL)
        {
            //printf(">>>>>>>>>>Eroare la deschiderea fisierului de intrare!<<<<<<<<<<");
            MPI_Abort(MPI_COMM_WORLD,-1);
            exit(-1);
        }
        //facem citirea datelor
        citireDateIntrare();
        //printf("-------------Afisare pentru citire date incheiata -------------\n");
        
        //printf("Master: lin,col: %d %d\n",lin,col);
    }

    //trimitem broadcast catre toate procesele cu numarul de linii si coloane al matricei initiale
    MPI_Bcast(&lin,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&col,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&valMax,1,MPI_INT,0,MPI_COMM_WORLD);

    //fiecare proces isi aloca matricea pentu fasia pe care o va prelucra
    if(rank!=nrProc-1)  //daca nu e vorba de ultimul proces
    {        
        nrLin=lin/nrProc+2;
        nrCol=col+2;
        //printf("nrProc,proc crt,lin,col: %d %d %d %d\n",nrProc,rank,lin,col);
        //printf("Calcul nrLin,nrCol pt >>rank %d<<: %d %d \n",rank,nrLin,nrCol);
        matrice=(int **) calloc(nrLin,sizeof(int*));
        for(i=0;i<nrLin;i++)
            matrice[i]=(int *) calloc(nrCol,sizeof(int));
    }
    else                //e vorba despre ultimul proces
    {
        nrLin=lin/nrProc+2+lin%nrProc;
        nrCol=col+2;
        //printf("nrProc,proc crt,lin,col: %d %d %d %d\n",nrProc,rank,lin,col);
        //printf("Calcul nrLin,nrCol pt >>rank %d<<: %d %d \n",rank,nrLin,nrCol);
        matrice=(int **) calloc(nrLin,sizeof(int*));
        for(i=0;i<nrLin;i++)
            matrice[i]=(int *) calloc(nrCol,sizeof(int));
    }

    if(rank==0) //daca e vorba de procesul root
    {        
        //isi copiaza in matrice liniile din matricea img care ii revin
        for(i=1;i<=nrLin-2;i++)
            for(j=1;j<=nrCol-2;j++)
                matrice[i][j]=img[i][j];    //-1 -1
        //debug
    /*    printf("--------------Matrice pentru root-------------\n");
        for(i=0;i<nrLin;i++)
        {
            for(j=0;j<nrCol;j++)
                printf("%d ",matrice[i][j]);
            printf("\n");
        }
        printf("--------------END Matrice pentru root-------------\n");
     */

        //trimite mesaje tuturor celorlalte procese cu liniile ce le revin acestora

        for(i=1;i<nrProc;i++)
        {
            //in functie de ce proces e calculam cate linii trebuiesc trimise
            if(i==nrProc-1) //e ultimul proces
                nrLinProc=nrLin+lin%nrProc;
            else
                nrLinProc=nrLin;

            //se trimit liniile corespunzatoare fasiei pe care o va prelucra procesul i
            for(j=i*(nrLin-2)+1;j<i*(nrLin-2)+(nrLinProc-2)+1;j++)
            {
                MPI_Send(img[j],nrCol,MPI_INT,i,DATE,MPI_COMM_WORLD);
                //debug
                /*printf("Trimis catre >proces %d< : ",i);
                for(k=0;k<nrCol;k++)
                    printf("%d ",img[j][k]);
                printf("\n");*/
            }
        }
        
    }

    if(rank>0)    //daca nu este vorba de procesul root
    {
        //fiecare proces primeste liniile care ii revin din matricea initiala
        for(i=1;i<=nrLin-2;i++)
            MPI_Recv(matrice[i],nrCol,MPI_INT,0,DATE,MPI_COMM_WORLD,&status);

        //debug
        /*printf("--------------Matrice pentru rank %d-------------\n",rank);
        for(i=0;i<nrLin;i++)
        {
            for(j=0;j<nrCol;j++)
                printf("%d ",matrice[i][j]);
            printf("\n");
        }
        printf("--------------END Matrice pentru rank %d-------------\n",rank);    */
    }

    //verificam care este tipul de ajustare pe crae il facem imaginii
    if(atoi(argv[1])==0)    //contrast
    {
        //fiecare proces face procesare informatiilor din fasia sa pentru a obtine portiunea modificata din imagine
        contrastImagine(rank,nrLin,nrCol,matrice,atoi(argv[3]),atoi(argv[4]));
    }
    else
        if(atoi(argv[1])==1)    //filtru
        {
            filtruImagine(rank,nrLin,nrCol,matrice,argv[3]);
        }
        else
            if(atoi(argv[1])==2) //entropie
            {
                entropieImagine(rank,nrLin,nrCol,matrice,atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
            }
            else    //nu e specificat corect
            {
                printf(">>>>>>>>>>Parametrul dat nu este corect! 0-contrast; 1-filtru; 2-entropie <<<<<<<<<<");
                MPI_Abort(MPI_COMM_WORLD,-1);
                exit(-1);
            }


    MPI_Barrier(MPI_COMM_WORLD);    //se asteapta sa termine toate de prelucrat

    if(rank>0)
    {
        //trimit raspunsul
        for(i=1;i<=nrLin-2;i++)
        {
            MPI_Send(matrice[i],nrCol,MPI_INT,0,SOL,MPI_COMM_WORLD);
            //debug
            //printf("\nproc <%d> a trimis linia [%d]\n",rank,i);
        }
        //debug
        /*printf("--------------Matricea modificata trimisa de proc %d-------------\n",rank);
        for(i=0;i<nrLin;i++)
        {
            for(j=0;j<nrCol;j++)
                printf("%d ",matrice[i][j]);
            printf("\n");
        }
        printf("--------------END Matricea modificata trimisa de proc %d-------------\n",rank);*/
    }

    //procesul root asteapta rezultatele si le utilizeaza pentru a-si forma noua matrice imagine
    if(rank==0)
    {
        int *proces;            //proces[i] va retine cate linii au fost primite de la procesul i
        int linPrimite=nrLin-2; //numarul de linii primite pana acum (momentan sunt cele ale procesului master(root))
        int pr;                 //va retine de la e proces e linia primita curent
        int linie;              //linia pe care trebuie copiata linia primita in noua imagine
        proces=(int *) calloc(nrProc,sizeof(int));

        //primeste informatii de la procese si in functie de acestea reinoieste imaginea
            //copiem matricea din procesul 0 in noua imagine
            for(i=1;i<=nrLin-2;i++)
                for(j=1;j<=nrCol-2;j++)
                    img[i][j]=matrice[i][j];

            //primeste o linie din matricea altui proces
            while(linPrimite<lin)
            {                
                MPI_Recv(matrice[1],nrCol,MPI_INT,MPI_ANY_SOURCE,SOL,MPI_COMM_WORLD,&status);
                
                linPrimite++;

                pr=status.MPI_SOURCE;
                proces[pr]++;    //se inoiste numarul de linii ce au fost primite de la procesul respectiv

                //aflam pe ce linie trebuie copiata linia primita in noua imagine
                linie=pr*(nrLin-2)+proces[pr];

                //debug
                //printf("!! S-a primit linia [[ %d ]] din procesul << %d >>\n",proces[pr],pr);

                //copiem linia respectiva in noua imagine
                for(j=1;j<=nrCol-2;j++)
                    img[linie][j]=matrice[1][j];
            }

        //daca este cazul in care trebuie obtinuta matricea reziduala / si calcularea entropiei
        if(atoi(argv[1])==2)
        {
            int minim=MAX_POSIBIL;
            int maxim=MIN_POSIBIL;
            double entropie=0;

            //aflam valoarea minima, respectiv maxima din matrice
            for(i=1;i<=lin;i++)
                for(j=1;j<=col;j++)
                    if(img[i][j]>maxim)
                        maxim=img[i][j];
                    else
                        if(img[i][j]<minim)
                            minim=img[i][j];

            //valorile ce apar vor fi retinute in vectSimb, nr lor de aparitii in vectAp, iar probabilitatea de aparitie se va calcula in vectProb
            int *vectSimb;
            vectSimb=(int *) calloc(maxim-minim+1,sizeof(int));

            int *vectAp;
            vectAp=(int *) calloc(maxim-minim+1,sizeof(int));

            float *vectProb;
            vectProb=(float *) calloc(maxim-minim+1,sizeof(float));

            //parcurgem matricea si aflam datele de care avem nevoie
            for(i=1;i<=lin;i++)
                for(j=1;j<=col;j++)
                {
                    vectSimb[img[i][j]-minim]=img[i][j];    //retinem ce valoare a aparut la pozitia -minim+valoare (pentru a nu )
                    vectAp[img[i][j]-minim]++;
                }
            for(i=0;i<maxim-minim+1;i++)
                if(vectAp[i]!=0)
                    vectProb[i]=(vectAp[i]+0.0)/(lin*col);

            //debug
            //printf("------------>SIMBOLURI, APARITII, PROBABILITATI<---------------\n");
            /*for(i=0;i<maxim-minim+1;i++)
                if(vectAp[i]!=0)
                {
                    printf("numar, nr aparitii, probabilitate: %d, %d, %f\n",vectSimb[i],vectAp[i],vectProb[i]);
                }
             */
            //printf("------------>END SIMBOLURI, APARITII, PROBABILITATI<---------------\n");
            //end debug

            for(i=0;i<maxim-minim+1;i++)
                if(vectAp[i]!=0)
                {
                    entropie=entropie-vectProb[i]*log2(vectProb[i]);
                }

            printf("Entropie: %f\n",entropie);

            //eliberam memoria alocata
            free(vectAp);
            free(vectProb);
            free(vectSimb);

            //inchidem fisierul de out si il redeschidem pentru a sterge ce mai fusese scris in el
            fclose(fout);
            fout=fopen(argv[6],"w");

            //afisam numarl de linii si coloane
            fprintf(fout,"%d %d\n",col,lin);
        }
        
        //debug
        //printf("----------------------NOUA MATRICE-----------------\n");
        //se face afisarea in fisier
        for(i=1;i<=lin;i++)
        {
            for(j=1;j<=col;j++)
                fprintf(fout,"%d\n",img[i][j]);
            //fprintf(fout,"%s","\n");
        }
        fprintf(fout,"%s","\n");
        //printf("---------------------- END NOUA MATRICE-----------------\n");

        //eliberam memoria alocata pentru a retine cate linii au fost primite de la fiecare proces
        free(proces);
    }

    //fiecare proces face dezalocare e memorie
    for(i=0;i<nrLin;i++)
        free(matrice[i]);
    free(matrice);

    //facem dezalocare de memorie pentru img
    if(rank==0)
    {
        for(i=0;i<lin+2;i++)
            free(img[i]);
        free(img);

        fclose(fin);
        fclose(fout);
    }

    //debug
    //printf("\nlin,col: %d %d\n",lin,col);

    MPI_Finalize();
}
