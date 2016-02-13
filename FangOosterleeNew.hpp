FangOost::FangOost(double xmin_, double xmax_, int k_, int h_){
  set_values(xmin_, xmax_, k_, h_);
}

void FangOost::set_values(double xmin_, double xmax_, int k_, int h_){
  xmin=xmin_;
  xmax=xmax_;
  h=h_;
  k=k_;
}
void FangOost::computeInv(auto&& fnInv){
  double xRange=xmax-xmin;
  du=M_PI/xRange;
  dx=xRange/(double)(h-1);
  double cp=2.0/xRange;
  f=std::vector<double> (k);
  #pragma omp parallel//multithread using openmp
  {
    #pragma omp for //multithread using openmp
    for(int j=0; j<k; j++){
      Complex u=Complex(0, du*j);
      f[j]=fnInv(u).multiply(u.multiply(-xmin).exp()).getReal()*cp;
    }
  }
  f[0]=.5*f[0];
}
std::vector<double> FangOost::computeConvolution(auto& vK){ //vk as defined in fang oosterlee
  std::vector<double> y(h);
  for(int i=0;  i<h; i++){
    y[i]=0;
    for(int j=0; j<k; j++){
      y[i]=y[i]+f[j]*vK(du*j, dx*i, xmin);//*cos(du*j*dx*i);
    }
  }
  return y;
}
void FangOost::computeExpectation(auto& callback){
  callback(f, du, dx, k, h, xmin);//accept f as parameter, will do something with f
}
