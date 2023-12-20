#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iomanip>

typedef double real;
typedef double realscal;
int precision = 14;
int n;
int notNullElements;
using namespace std;

void inputSize(string path){
   ifstream file(path);
   if(file.is_open()){
      file >> n;
      file.close();
   }
   else{
      cerr << "Can't open the input file" << path << endl;
   }
}

void inputNotNullElements(string path){
   ifstream file(path);
   if(file.is_open()){
      while(file >> notNullElements){};
      notNullElements--;
      file.close();
   }
   else{
      cerr << "Can't open the input file" << path << endl;
   }
}

struct sole{
   int shift = 0;
   real *diagnal = new real[n];
   real *f = new real[n];
   real *y = new real[n];
   real *x = y;
   int *ia = new int[n];
   real *au = new real[notNullElements];
   real *al = new real[notNullElements];

   void inputDiagnal(string path);
   void inputIa(string path);
   void inputData(string path, real *arr);
   void inputVector(string path, real *arr);
   void calculateL(int i, int j);
   void calculateU(int i, int j);
   void calculateY();
   void calculateX();
   void LUdecomposition();
   void printAllData();
   void printLU();
   void outputX(string path);
   void program1();
   void program2();
   void Gauss();
   void program3();
   void prm(real m[][10]);
};

void sole::inputDiagnal(string path){
   ifstream file(path);
   if(file.is_open()){
      for(int i = 0; i < n; i++){
         file >> diagnal[i];
      }
      file.close();
   }else{
      cerr << "Can't open the input file" << path << endl;
   }
}

void sole::inputIa(string path){
   ifstream file(path);
   if(file.is_open()){
      for(int i = 0; i<n+1; i++){
         file >> ia[i];
      }
      file.close();
   }
   else{
      cerr << "Can't open the input file" << path << endl;
   }
}

void sole::inputData(string path, real *arr){
   ifstream file(path);
   if(file.is_open()){
      for(int i = 0; i < notNullElements; i++){
         file >> arr[i];
      }
      file.close();
   }
   else{
      cerr << "Can't open the input file" << path << endl;
   }
}

void sole::inputVector(string path, real *arr){
   ifstream file(path);
   if(file.is_open()){
      for(int i = 0; i < n; i++){
         file >> arr[i];
      }
      file.close();
   } else{
      cerr << "Can't open the input file" << path << endl;
   }
}


void sole::printAllData(){
   for(int i = 0; i<n+1; i++){
      cout << ia[i] << endl;
   }cout << "\n\n\n";

   for(int i = 0; i<notNullElements; i++){
      cout << al[i] << endl;
   }cout << "\n\n\n";

   for(int i = 0; i<notNullElements; i++){
      cout << au[i] << endl;
   }cout << "\n\n\n";

   for(int i = 0; i<n; i++){
      cout << diagnal[i] << endl;
   }cout << "\n\n\n";
}

void sole::printLU(){
   cout << "U matrix:" << endl;
   for(int i = 0; i<n; i++){
      for(int j = 0; j<n; j++){
         if(i==j){
            cout << diagnal[i] << setw(16) << "\t";
         } else{
            if(i>j || ia[j+1]-ia[j] == 0 || j-(ia[j+1]-ia[j]) > i){
               cout << 0 << setw(16) << "\t";
            } else{
               cout << au[ia[j+1]-1-(j-i)] << setw(16) << "\t";
            }
         }
      }
      cout << endl;
   }cout << "\n\n\n";

   cout << "L matrix:" << endl;
   for(int i = 0; i<n; i++){
      for(int j = 0; j<n; j++){
         if(i==j){
            cout << 1 << setw(16) << "\t";
         } else{
            if(i<j || ia[i+1]-ia[i] == 0 || i-(ia[i+1]-ia[i]) > j){
               cout << 0 << setw(16) << "\t";
            } else{
               cout << al[ia[i+1]-1-(i-j)] << setw(16) << "\t";
            }
         }
      }
      cout << endl;
   }cout << "\n\n\n";
}

void sole::calculateU(int i, int j){
   realscal sum = 0;
   real oldElem = 0, L = 0, U = 0;
   if(i==j){
      oldElem = diagnal[i];
   } else{
      if(j-shift > i){
         oldElem = 0;
      } else{
         oldElem = au[ia[j+1]-1-(j-i)];
      }
   }

   for(int k = 0; k<=i-1 && k!=j; k++){
      if(i==k){
         L = 1;
      }
      else{
         if(ia[i+1]-ia[i] == 0 || i-(ia[i+1]-ia[i]) > k){
            L = 0;
         } else{
            L = al[ia[i+1]-1-(i-k)];
         }
      }

      if(j-shift > k){
         U = 0;
      }
      else{
         U = au[ia[j+1]-1-(j-k)];
      }
   sum = sum + L * U;
   }

   if(i == j){ 
      sum = oldElem - sum;
      diagnal[i] = sum;
   }
   else{
      sum = oldElem - sum;
      au[ia[j+1]-1-(j-i)] = sum;
   }
}

void sole::calculateL(int i, int j){
   realscal sum = 0;
   real oldElem = 0, L = 0, U = 0;
   if(i-shift > j){
      oldElem = 0;
   } else{
      oldElem = al[ia[i+1]-1-(i-j)];
   }

   for(int k = 0; k<=j-1 && k!=j; k++){
      if(i==k){
         L = 1;
      } else{
         if(ia[i+1]-ia[i] == 0 || i-(ia[i+1]-ia[i]) > k){
            L = 0;
         } else{
            L = al[ia[i+1]-1-(i-k)];
         }
      }

      if(ia[j+1]-ia[j] == 0 || j-(ia[j+1]-ia[j]) > k){
         U = 0;
      } else{
         U = au[ia[j+1]-1-(j-k)];
      }
      sum = sum+L*U;
   }
   sum = (oldElem-sum)/diagnal[j];
   al[ia[i+1]-1-(i-j)] = sum;
}

void sole::LUdecomposition(){
   for(int h = 0; h<n; h++){
      for(int j = h; j<=n-1; j++){
         shift = ia[j+1]-ia[j];

         if(h==j){
            calculateU(h, j);
         } else{
            if(shift != 0 && j-shift <= h){ 
               calculateU(h, j);
            }
         }
      }

      for(int i = h+1; i<n; i++){
         shift = ia[i+1]-ia[i];

         if(shift != 0 && i-shift <= h){
            calculateL(i, h);
         }
      }
   }
}

void sole::calculateY(){
   realscal sum;
   real L;
   for(int i = 0; i<n; i++){
      sum = 0;
      for(int k = 0; k<=i-1; k++){
         if(i==k){
            L = 1;
         } else{
            if(ia[i+1]-ia[i] == 0 || i-(ia[i+1]-ia[i]) > k){
               L = 0;
            } else{
               L = al[ia[i+1]-1-(i-k)];
            }
         }
         sum = sum + L*y[k];
      }
      y[i] = f[i]-sum;
   }
}

void sole::calculateX(){
   realscal sum;
   real U;
   for(int i = n-1; i!= -1; i--){
      sum = 0;
      for(int k = n-1; k>i; k--){
         if(i==k){
            U = diagnal[i];
         } else{
            if(ia[k+1]-ia[k] == 0 || k-(ia[k+1]-ia[k]) > i){
               U = 0;
            } else{
               U = au[ia[k+1]-1-(k-i)];
            }
         }

         sum = sum + U*x[k];
      }
      x[i] = (y[i]-sum)/diagnal[i];
   }
}

void sole::outputX(string path){
   ofstream file(path);
   if(file.is_open()){
      for(int i = 0; i<n-1; i++){
         file << setprecision(precision) << x[i] << endl;
      }file << setprecision(precision) << x[n-1];
      file.close();
   } else{
      cerr << "Can't open the output file" << path << endl;
   }

}

void sole::program1(){   //LU декомпозиция и вычисление векторов поочерёдно
   inputDiagnal("diagnal.txt");
   inputIa("ia.txt");
   inputData("al.txt", al);
   inputData("au.txt", au);
   inputVector("f.txt", f);
   LUdecomposition();
   calculateY();
   calculateX();
   outputX("output.txt");
}


void generateHilbert(){
   ofstream iaFile("ia.txt");
   ofstream alFile("al.txt");
   ofstream auFile("au.txt");
   ofstream diagnalFile("diagnal.txt");
   realscal sum = 0;
   int index = 0;
   diagnalFile << 1;
   iaFile << 1 << endl << 1;
   if(n>1){ index = 1; }
   for(int i = 1; i<n; i++){
      diagnalFile << endl << pow(2*i+1,-1);
      for(int j = 0; j<i; j++){
         sum = pow(i+j+1, -1);
         alFile << sum << endl;
         auFile << sum << endl;
         index++;
      }
      iaFile << endl << index;
   }
   notNullElements = index-1;
   iaFile.close();
   alFile.close();
   auFile.close();
   diagnalFile.close();
}

void sole::program2(){   //LU декомпозиция и вычисление векторов поочерёдно сгенерированной матрицы Гильберта
   inputDiagnal("diagnal.txt");
   inputIa("ia.txt");
   inputData("al.txt", al);
   inputData("au.txt", au);
   inputVector("f.txt", f);
   printAllData();
   LUdecomposition();
   printLU();
   printAllData();
   calculateY();
   calculateX();
   outputX("output.txt");
}

void sole::Gauss(){
   real matrix[10][11]{};

   for(int i = 0; i<n; i++){
      for(int j = 0; j<n; j++){
         if(i==j){
            matrix[i][j] = diagnal[i];
         } else{
            if(i<j){
               if(ia[j+1]-ia[j] != 0 && j-(ia[j+1]-ia[j]) <= i){
                  matrix[i][j] = au[ia[j+1]-1-(j-i)];
               }
            } else{
               if(ia[i+1]-ia[i] != 0 && i-(ia[i+1]-ia[i]) <= j){
                  matrix[i][j] = al[ia[i+1]-1-(i-j)];
               }
            }
         }
      }
      matrix[i][n] = f[i];
   }

   for(int k = 0; k<n; k++){
      for(int i = 0; i<=n; i++){
         matrix[k][i] = matrix[k][i] / diagnal[k];
      }
      for(int i = k+1; i<n; i++){
         real K = matrix[i][k] / matrix[k][k];
         for(int j = 0; j<=n; j++){
            matrix[i][j] = matrix[i][j] - matrix[k][j] * K;
         }
      }
   }

   for(int k = n-1; k>=0; k--){
      real K = 1/matrix[k][k];
      for(int j = n; j>=k; j--){
         matrix[k][j] *= K;
      }
   }

   for(int k = n-1; k>=0; k--){
      real buf = 0;
      for(int j = n-1; k!= n-1 && j>k; j--){
         buf += matrix[k][j]*x[j];
      }
      x[k] = matrix[k][n] - buf;
   }
}

void sole::program3(){   //Вычисление вектора х по методу Гаусса из матрицы размерности 10
   inputDiagnal("diagnal.txt");
   inputIa("ia.txt");
   inputData("al.txt", al);
   inputData("au.txt", au);
   inputVector("f.txt", f);
   Gauss();
   outputX("output.txt");
}

void main()
{
   inputSize("size.txt");
   //generateHilbert();   
   inputNotNullElements("ia.txt");
   sole data;
   data.program3();
}