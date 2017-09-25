#include <stdio.h>
#include <getopt.h> // *GNU* Para o getopt_long()
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(int argc, char *argv[])
{
FILE *vasprun,*imag,*real,*dost,*aux;
int i,j,nedos,plot,chi,nomega;
char str[200],ch;
float pi,hbar,e,e0,number,energy,c;

hbar=4.13566743*pow(10,-15);  /*eV.s*/
/*hbar=1.05459*pow(10,-34);  /*J.s*/ 
pi=3.14159;
e=1.60219*pow(10,-19); /*coulombs*/
e0=8.85*pow(10,-12); /*permissividade vacuo, epsilon zero C^2/(Nm^2)*/
c=299792458; /*velocidade da luz em m/s*/


/*-----------arguments---------------*/

	int opt=0,wave=0;  

//Specifying the expected options

static struct option long_options[] = {
        {"twodim",      no_argument,       0,  't' },
        {"wavelenght", no_argument,       0,  'w' },
        {"equations", no_argument,       0,  'e' },
        {"help",    no_argument, 0,  'h' },
        {0,           0,                 0,  0   }
    };


int long_index =0;
    while ((opt = getopt_long(argc, argv,"tweh",long_options, &long_index )) != -1) {
        switch (opt) {
                        case 't' : // Usuario
                                printf("For two dimensional systems! not implemented yet\n");
                                exit(0);
                                break;
                        case 'w' : // Arquivo
                                wave=1;printf("\nREMEMBER TO ADJUST THE SCALE OF X-AXE!\n\n");
                                break;
                        case 'e' : // Arquivo
                                printf("\nDefinitions\n\n");
                                printf("\nRe[e]: Real part of dielectric function");
                                printf("\nIm[e]: Imaginaty part of dielectric function");
                                printf("\nomega: Frequency of the light");
                                printf("\nE: Energy of the light in eV");
                                printf("\nc: light velocity in m/s");
                                printf("\nhbar: h/2pi in ev.s");
                                printf("\ne0: Permitivity of vacuum in C^2/(Nm^2)");
                                printf("\ne: Charge of electron in Coulomb\n");

                                printf("\nAbsorption coefficient: alpha = 2*omega*k/c = 2*E*k/c/hbar");
                                printf("\nRefractive index: n = sqrt((sqrt(Re[e]^2+Im[e]^2)+Re[e])/2)");
                                printf("\nExtinction coefficient: k = sqrt((sqrt(Re[e]^2+Im[e]^2)-Re[e])/2)");
                                printf("\nSOLID STATE PHYSICS PART II Optical Properties of Solids M. S. Dresselhaus\n");
                                printf("\nOptical conductivity: sigma = E*Im[e]*4*e0/e");
                                printf("\nNew Journal of Physics 16 (2014) 105007\n");
                                printf("\nReflectance: R = r=((n-1)^2+k^2)/((n+1)^2+k^2)");
                                printf("\n2013 J. Semicond. 34 032002\n\n");exit(0);
                                
                        case 'h' : // Ajuda
                                printf("\nUse opticsVASP without any options to get default values or:\n");
                                printf("\n-t For two dimensional systems, including a vacuum part in z direction!\n");
                                printf("-w x-axe in wavenumber (nm) insted eV.\n");
                                printf("-e to see the equations used.\n");
                                printf("-h print this help.\n\n");
                                exit(0);
                        default : // Qualquer parametro nao tratado
                                printf("Incorrect parameter, please use option -h for help!\n");
                                exit(0);
        }
    }


printf("\n\n");
printf("Program to plot optical properties from VASP results\n");
printf("by F. Matusalem - 07/2017 - filipematus@gmail.com\n");
printf("Version 1.0 07/2017\n\n");

printf("use -h option for help\n\n");

printf("Absorption coefficient (0)\n"); /*SOLID STATE PHYSICS PART II Optical Properties of Solids M. S. Dresselhaus*/
printf("      Refractive index (1)\n");
printf("Extinction coefficient (2)\n");
printf("  Optical conductivity (3)\n");  /*New Journal of Physics 16 (2014) 105007*/
printf("   Imag. of dielectric (4)\n");
printf("    Real of dielectric (5)\n");
printf("           Reflectance (6)\n");  /*2013 J. Semicond. 34 032002*/

printf("\nChoose the optical property to plot: \n");
scanf("%d",&plot);



vasprun = fopen("vasprun.xml","r"); /* Arquivo ASCII, para leitura */
if(!vasprun)
{
printf( "Erro na abertura do arquivo vasprun.xml\n");
exit(0);
}


imag = fopen("imag-dielectric.dat","w"); /* Arquivo ASCII, para escrita */
if(!imag)
{
printf( "Error creating imag-dielectric.dat file\n");
exit(0);
}

real = fopen("real-dielectric.dat","w"); /* Arquivo ASCII, para escrita */
if(!real)
{
printf( "Error creating real-dielectric.dat file\n");
exit(0);
}

i=0;
while (fscanf(vasprun,"%s",str) != EOF){
if(strcmp(str,"name=\"ALGO\">")==0)i++;                      /*verifica se o arquivo vasprun tem name=\"ALGO\">*/
}
rewind(vasprun);

chi=0;
if(i!=0){
do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra algo*/
while(strcmp(str,"name=\"ALGO\">")!=0);
fscanf(vasprun,"%s",str);      /*lê o algo*/
if(strcmp(str,"chi</i>")==0)chi=1;
}


do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra NEDOS*/
while(strcmp(str,"name=\"NEDOS\">")!=0);
fscanf(vasprun,"%d",&nedos);      /*lê o numero de NEDOS*/


do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra NEDOS*/
while(strcmp(str,"name=\"NOMEGA\">")!=0);
fscanf(vasprun,"%d",&nomega);      /*lê o numero de name="NOMEGA">*/



float img[nedos][7],re[nedos][7],k,l,n,r;



if(chi==0){


i=0;
while (fscanf(vasprun,"%s",str) != EOF){
if(strcmp(str,"<dielectricfunction>")==0)i++;                      /*verifica se o arquivo vasprun tem a função dieletrica*/
}

if(i == 0){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");


printf("\n\nThe vasprun.xml file do not have the dielectric function. Rerun the calculation using LOPTICS = T!!\n");
printf("\nIf you are already using LOPTICS = T, try to rerun using ALGO = nothing ; LPEAD = .TRUE.!! bye! bye! \n\n");
exit(0);} 
rewind(vasprun);


do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<dielectricfunction>")!=0);

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<imag>")!=0);

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<set>")!=0);

do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');


for(i=0;i<nedos;i++){fscanf(vasprun,"%s",str);
for(j=0;j<7;j++)fscanf(vasprun,"%f",&img[i][j]);
do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(i=0;i<nedos;i++){for(j=0;j<7;j++)fprintf(imag,"%f ",img[i][j]);fprintf(imag,"\n");}

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<real>")!=0);

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<set>")!=0);


do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');

for(i=0;i<nedos;i++){fscanf(vasprun,"%s",str);
for(j=0;j<7;j++)fscanf(vasprun,"%f",&re[i][j]);
do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(i=0;i<nedos;i++){for(j=0;j<7;j++)fprintf(real,"%f ",re[i][j]);fprintf(real,"\n");}

fclose(vasprun);
fclose(imag);
fclose(real);
}

else {

nedos=nomega;

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra MACROSCOPIC*/
while(strcmp(str,"MACROSCOPIC")!=0);

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<imag>")!=0);

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<set>")!=0);

do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');



for(i=0;i<nedos;i++){fscanf(vasprun,"%s",str);
for(j=0;j<7;j++)fscanf(vasprun,"%f",&img[i][j]);
do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(i=0;i<nedos;i++){for(j=0;j<7;j++)fprintf(imag,"%f ",img[i][j]);fprintf(imag,"\n");}

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<real>")!=0);

do
fscanf(vasprun,"%s",str);                                      /*posiciona o ponteiro após a palavra IMAGINARY*/
while(strcmp(str,"<set>")!=0);


do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');

for(i=0;i<nedos;i++){fscanf(vasprun,"%s",str);
for(j=0;j<7;j++)fscanf(vasprun,"%f",&re[i][j]);
do
ch = getc(vasprun);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(i=0;i<nedos;i++){for(j=0;j<7;j++)fprintf(real,"%f ",re[i][j]);fprintf(real,"\n");}

fclose(vasprun);
fclose(imag);
fclose(real);
}









/*-------------------------------------------------------------------------------------------------------*/

switch ( plot )
  {
     case 0 :
     
aux = fopen("abs-coefficient.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo abs-coefficient.dat\n");
exit(0);
}

if(wave==0){
for(i=0;i<nedos;i++){
fprintf(aux,"%f ",re[i][0]);
for(j=1;j<7;j++){
k=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))-re[i][j])/2);

fprintf(aux,"%f ",2*re[i][0]*k/c/hbar/10000000);}fprintf(aux,"\n");
}
}

else {
float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++){
k=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))-re[i][j])/2);

fprintf(aux,"%f ",2*re[i][0]*k/c/hbar/10000000);}fprintf(aux,"\n");
}



}

fclose(aux);



dost = fopen("abs-coefficient.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo abs-coefficient.bfile\n");
exit(0);
}

fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
fprintf(dost,"read block \"abs-coefficient.dat\"\n");
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Absorption coefficient (10\\S7 \\N cm\\S-1\\N)\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param abs-coefficient.bfile -saveall abs-coefficient.agr -hdevice EPS");
break;


     case 1 :
     
aux = fopen("refractiveindex.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo refractiveindex.dat\n");
exit(0);
}

if(wave==0){
for(i=0;i<nedos;i++){
fprintf(aux,"%f ",re[i][0]);
for(j=1;j<7;j++){
n=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))+re[i][j])/2);

fprintf(aux,"%f ",n);}fprintf(aux,"\n");
}}


else {
float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++){
n=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))+re[i][j])/2);

fprintf(aux,"%f ",n);}fprintf(aux,"\n");
}}



fclose(aux);



dost = fopen("refractiveindex.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo refractiveindex.bfile\n");
exit(0);
}

fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
fprintf(dost,"read block \"refractiveindex.dat\"\n");
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Refractive Index\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param refractiveindex.bfile -saveall refractiveindex.agr -hdevice EPS");
break;




     case 2 :
     
aux = fopen("extinctioncoefficient.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo extinctioncoefficient.dat\n");
exit(0);
}

if(wave==0){
for(i=0;i<nedos;i++){
fprintf(aux,"%f ",re[i][0]);
for(j=1;j<7;j++){
l=img[i][j]/(pow(re[i][j],2)+pow(img[i][j],2));
fprintf(aux,"%f ",l);}fprintf(aux,"\n");
}}

else {
float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++){
l=img[i][j]/(pow(re[i][j],2)+pow(img[i][j],2));
fprintf(aux,"%f ",l);}fprintf(aux,"\n");
}}









fclose(aux);



dost = fopen("extinctioncoefficient.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo extinctioncoefficient.bfile\n");
exit(0);
}

fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
fprintf(dost,"read block \"extinctioncoefficient.dat\"\n");
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Energy Loss Spectrum\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param extinctioncoefficient.bfile -saveall extinctioncoefficient.agr -hdevice EPS");
break;


     case 3 :
     
aux = fopen("opticalconductivity.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo opticalconductivity.dat\n");
exit(0);
}
if(wave==0){
for(i=0;i<nedos;i++){
fprintf(aux,"%f ",re[i][0]);
for(j=1;j<7;j++)fprintf(aux,"%f ",re[i][0]*img[i][j]*4*e0/e);
fprintf(aux,"\n");
}}


else {
float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++)fprintf(aux,"%f ",re[i][0]*img[i][j]*4*e0/e);
fprintf(aux,"\n");
}}

fclose(aux);



dost = fopen("opticalconductivity.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo opticalconductivity.bfile\n");
exit(0);
}

fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
fprintf(dost,"read block \"opticalconductivity.dat\"\n");
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Optical conductivity Re \\xs\\f{}\\s2D(\\xw)\\N/\\xs\\s0\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param opticalconductivity.bfile -saveall opticalconductivity.agr -hdevice EPS");
break;


     case 4 :
     
dost = fopen("imag-dielectric.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo imag-dielectric.bfile\n");
exit(0);
}

if(wave!=0){

aux = fopen("imag-dielectric-wave.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo imag-dielectric-wave.dat\n");
exit(0);
}

float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++){fprintf(aux,"%f ",img[i][j]);fprintf(aux,"\n");}
}

fclose(aux);}



fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
if(wave==0)fprintf(dost,"read block \"imag-dielectric.dat\"\n");
else {fprintf(dost,"read block \"imag-dielectric-wave.dat\"\n");}
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Im \\xe[\\xw]\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param imag-dielectric.bfile -saveall imag-dielectric.agr -hdevice EPS");
break;

     case 5 :
     
dost = fopen("real-dielectric.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo real-dielectric.bfile\n");
exit(0);
}

if(wave!=0){
aux = fopen("real-dielectric-wave.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo real-dielectric-wave.dat\n");
exit(0);
}

float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++){fprintf(aux,"%f ",re[i][j]);fprintf(aux,"\n");}
}
fclose(aux);}




fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
if(wave==0)fprintf(dost,"read block \"real-dielectric.dat\"\n");
else {fprintf(dost,"read block \"real-dielectric-wave.dat\"\n");}
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Re \\xe[\\xw]\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param real-dielectric.bfile -saveall real-dielectric.agr -hdevice EPS");
break;

     case 6 :
     
aux = fopen("reflectance.dat","w"); /* Arquivo ASCII */
if(!aux)
{printf( "Erro na criacao do arquivo reflectance.dat\n");
exit(0);
}

if(wave==0){
for(i=0;i<nedos;i++){
fprintf(aux,"%f ",re[i][0]);
for(j=1;j<7;j++){
n=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))+re[i][j])/2);
k=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))-re[i][j])/2);
r=(pow(n-1,2)+pow(k,2))/(pow(n+1,2)+pow(k,2));

fprintf(aux,"%f ",r);}fprintf(aux,"\n");
}}


else {
float lambda;
for(i=1;i<nedos;i++){
lambda=2*pi*hbar*c/re[i][0];
fprintf(aux,"%f ",1000000000*lambda);
for(j=1;j<7;j++){
n=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))+re[i][j])/2);
k=sqrt((sqrt(pow(re[i][j],2)+pow(img[i][j],2))-re[i][j])/2);
r=(pow(n-1,2)+pow(k,2))/(pow(n+1,2)+pow(k,2));

fprintf(aux,"%f ",r);}fprintf(aux,"\n");
}}



fclose(aux);



dost = fopen("reflectance.bfile","w"); /* Arquivo ASCII */
if(!dost)
{printf( "Erro na criacao do arquivo reflectance.bfile\n");
exit(0);
}

fprintf(dost,"title font 14 \n");
fprintf(dost,"subtitle font 14  \n");
fprintf(dost,"legend font 14  \n");
fprintf(dost,"read block \"reflectance.dat\"\n");
fprintf(dost,"block xy \"1:2\" \n");
fprintf(dost,"s0 symbol 0 \n");
fprintf(dost,"s0 line type 1 \n");
fprintf(dost,"s0 line linestyle 1 \n");
fprintf(dost,"s0 line linewidth 1 \n");
fprintf(dost,"s0 line color 1 \n");
fprintf(dost,"s0 comment \"xx\" \n");
fprintf(dost,"s0 legend  \"xx\" \n");
fprintf(dost,"block xy \"1:3\" \n");
fprintf(dost,"s1 symbol 0 \n");
fprintf(dost,"s1 line type 1 \n");
fprintf(dost,"s1 line linestyle 1 \n");
fprintf(dost,"s1 line linewidth 1 \n");
fprintf(dost,"s1 line color 1 \n");
fprintf(dost,"s1 comment \"yy\" \n");
fprintf(dost,"s1 legend  \"yy\" \n");
fprintf(dost,"block xy \"1:4\" \n");
fprintf(dost,"s2 symbol 0 \n");
fprintf(dost,"s2 line type 1 \n");
fprintf(dost,"s2 line linestyle 1 \n");
fprintf(dost,"s2 line linewidth 1 \n");
fprintf(dost,"s2 line color 2 \n");
fprintf(dost,"s2 comment \"zz\" \n");
fprintf(dost,"s2 legend  \"zz\" \n");
fprintf(dost,"block xy \"1:5\" \n");
fprintf(dost,"s3 symbol 0 \n");
fprintf(dost,"s3 line type 1 \n");
fprintf(dost,"s3 line linestyle 1 \n");
fprintf(dost,"s3 line linewidth 1 \n");
fprintf(dost,"s3 line color 2 \n");
fprintf(dost,"s3 comment \"xy\" \n");
fprintf(dost,"s3 legend  \"xy\" \n");
fprintf(dost,"block xy \"1:6\" \n");
fprintf(dost,"s4 symbol 0 \n");
fprintf(dost,"s4 line type 1 \n");
fprintf(dost,"s4 line linestyle 1 \n");
fprintf(dost,"s4 line linewidth 1 \n");
fprintf(dost,"s4 line color 3 \n");
fprintf(dost,"s4 comment \"yz\" \n");
fprintf(dost,"s4 legend  \"yz\" \n");
fprintf(dost,"block xy \"1:7\" \n");
fprintf(dost,"s5 symbol 0 \n");
fprintf(dost,"s5 line type 1 \n");
fprintf(dost,"s5 line linestyle 1 \n");
fprintf(dost,"s5 line linewidth 1 \n");
fprintf(dost,"s5 line color 3 \n");
fprintf(dost,"s5 comment \"zx\" \n");
fprintf(dost,"s5 legend  \"zx\" \n");
fprintf(dost,"yaxis  tick on  \n");
fprintf(dost,"yaxis  tick major 4  \n");
fprintf(dost,"yaxis  tick minor ticks 1  \n");
fprintf(dost,"yaxis  tick major size 1.000000 \n");
fprintf(dost,"yaxis  ticklabel char size 1.200000 \n");
fprintf(dost,"yaxis  label \"Reflectance\" \n");
fprintf(dost,"yaxis  label font 0 \n");
fprintf(dost,"yaxis  label char size 1.200000 \n");
fprintf(dost,"xaxis  tick on  \n");
fprintf(dost,"xaxis  tick major 2  \n");
fprintf(dost,"xaxis  tick minor ticks 1  \n");
fprintf(dost,"xaxis  tick major size 1.000000 \n");
fprintf(dost,"xaxis  ticklabel char size 1.200000 \n");
if(wave==0)fprintf(dost,"xaxis  label \"Photon energy \\1h\\v{0.65}\\h{-0.5}\\z{0.6}_\\v{}\\h{}\\z{}\\xw \\f{}(eV)\" \n");
else {fprintf(dost,"xaxis  label \"Wavelenght (nm)\" \n");}
fprintf(dost,"xaxis  label font 0 \n");
fprintf(dost,"xaxis  label char size 1.200000 \n");
fclose(dost);

/*plota xmgrace*/
system("gracebat -param reflectance.bfile -saveall reflectance.agr -hdevice EPS");
break;


default :
       printf ("Valor invalido!\n");
  }


printf("\n\n"); 
printf("---------------------------------------------------------------------------------------\n"); 
printf("by F. Matusalem - 06/2017 version 2.0\n");
printf("Optics properties for VASP.5.3 results\n\n");
printf("---------------------------------------------------------------------------------------\n");   
}
