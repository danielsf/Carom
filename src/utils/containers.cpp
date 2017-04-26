#include "containers.h"

template <typename T>
asymm_array_2d<T>::asymm_array_2d(){
    name=NULL;
    where_am_i=NULL;
    data=NULL;
    rows=0;
    row_room=0;
}

template <typename T>
asymm_array_2d<T>::~asymm_array_2d(){
    int i;

    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name_null();
        }
        catch(int iex){
            printf("in asymm 2d destructor\n");
            die(0);
        }

        try{
            data[i].assert_where_null();
        }
        catch(int iex){
            printf("in asymm 2d destructor\n");
            die(0);
        }
    }

    //printf("calling 2d destructor on %s\n",name);

    if(data!=NULL){
        delete [] data;
    }

    if(name!=NULL)delete [] name;

    if(where_am_i!=NULL)delete [] where_am_i;


}

template <typename T>
void asymm_array_2d<T>::die(int ir) const{
    printf("\nWARNING asymm 2d array\n");

    if(name!=NULL)printf("in 2d_array %s\n",name);
    if(where_am_i!=NULL)printf("in routine %s\n",where_am_i);

    printf("asked for %d\n",ir);
    printf("but dimensions are %d\n",rows);
    printf("row_room %d\n",row_room);

    if(data==NULL){
        printf("data is null\n");
    }

    int ifail=1;

    throw ifail;
}

template <typename T>
void asymm_array_2d<T>::zero(){
    int i;
    for(i=0;i<row_room;i++)data[i].zero();
}

template <typename T>
void asymm_array_2d<T>::set_where(char *word) const {

    int i,ct=0;

    for(i=0;word[i]!=0;i++)ct++;
    ct++;

    if(where_am_i!=NULL)delete [] where_am_i;
    where_am_i=new char[ct];

    for(i=0;i<ct && word[i]!=0;i++)where_am_i[i]=word[i];
    where_am_i[ct-1]=0;

    for(i=0;i<row_room;i++){
        try{
            data[i].assert_where(where_am_i);
        }
        catch(int iex){
            printf("in asymm 2d set where\n");
            die(0);
        }
    }

}

template <typename T>
void asymm_array_2d<T>::set_name(char *word){
    int i,ct=0;

    for(i=0;word[i]!=0;i++)ct++;
    ct++;

    if(name!=NULL){
        delete [] name;
    }

    name=new char[ct];
    for(i=0;i<ct && word[i]!=0;i++)name[i]=word[i];
    name[ct-1]=0;

    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name(name);
        }
        catch(int iex){
            printf("in asymm 2d set name\n");
            die(0);
        }
    }

}

template <typename T>
void asymm_array_2d<T>::add_row(const array_1d<T> &in){

    if(data==NULL){
        row_room=2;
        rows=0;
        data=new array_1d<T>[row_room];
    }

    array_1d<T> *buffer;
    int i,j;

    if(rows==row_room){
        buffer=new array_1d<T>[rows];
        for(i=0;i<rows;i++){
            buffer[i].set_dim(data[i].get_dim());
            for(j=0;j<data[i].get_dim();j++){
                buffer[i].set(j,data[i].get_data(j));
            }
        }
        delete [] data;

        i=row_room/2;
        if(i<100)i=100;

        row_room*=2;
        data=new array_1d<T>[row_room];

        for(i=0;i<rows;i++){
            data[i].set_dim(buffer[i].get_dim());
            for(j=0;j<buffer[i].get_dim();j++){
                data[i].set(j,buffer[i].get_data(j));
            }
        }
        delete [] buffer;
    }

    data[rows].set_dim(in.get_dim());
    for(i=0;i<in.get_dim();i++){
        data[rows].set(i,in.get_data(i));
    }
    rows++;

    for(i=0;i<row_room;i++){
            try{
                data[i].assert_name(name);
            }
            catch(int iex){
                printf("in asymm 2d add row (asserting name)\n");
                die(0);
            }

            try{
                data[i].assert_where(where_am_i);
            }
            catch(int iex){
                printf("in asymm 2d add row (asserting where)\n");
                die(0);
            }
    }

}

template <typename T>
void asymm_array_2d<T>::set_row(int ir, const array_1d<T> &vv){

    int i;
    array_1d<T> empty;

    while(rows<=ir){
        add_row(empty);
    }

    data[ir].reset();
    for(i=0;i<vv.get_dim();i++){
        data[ir].set(i,vv.get_data(i));
    }


}

template <typename T>
void asymm_array_2d<T>::remove_row(int dex){

    if(dex<0 || dex>=rows){
        printf("WARNING asking to remove %d from asymm\n",dex);
        die(dex);
    }

    int i,j;
    for(i=dex;i<rows-1;i++){
        data[i].set_dim(data[i+1].get_dim());
        for(j=0;j<data[i+1].get_dim();j++){
            data[i].set(j,data[i+1].get_data(j));
        }
    }
    data[rows-1].reset();
    rows--;


}

template <typename T>
void asymm_array_2d<T>::reset_preserving_room(){
    int i;
    for(i=0;i<rows;i++){
        data[i].reset_preserving_room();
    }
    rows=0;
}

template <typename T>
void asymm_array_2d<T>::reset(){

    //printf("resetting %s\n",name);

    //set_name("resetting");

    int i;

    if(data==NULL && (rows>0 || row_room>0)){
        printf("resetting but data is null and something is wrong\n");
        die(-1);
    }

    if(row_room==0 && data!=NULL){
        die(-1);
    }

    if(row_room<rows){
        die(-2);
    }

    if(data!=NULL){
        delete [] data;
        data=NULL;
        row_room=0;
        rows=0;

    }



}

template <typename T>
void asymm_array_2d<T>::replace_row(int dex, array_1d<T> &pt){
    if(dex<0 || dex>=rows){
        printf("WARNING trying to replace row %d in asymm, but only have %d\n",
        dex,rows);

        die(dex);
    }

    data[dex].reset();
    int i;
    for(i=0;i<pt.get_dim();i++){
        data[dex].add(pt.get_data(i));
    }

}

template class asymm_array_2d<int>;
template class asymm_array_2d<double>;
template class array_2d<int>;
template class array_2d<double>;
template class array_1d<int>;
template class array_1d<double>;


template <typename T>
void merge_sort(array_1d<T> &in, array_1d<int> &dexes,
                int start, int end){

    if(end>=dexes.get_dim()){
        printf("WARNING in merge_sort end is %d but dexes has %d\n",
        end,dexes.get_dim());

        exit(1);
    }

    if(in.get_dim()==0) return;

    in.set_where("merge_sort");

    T nn;
    int i1,i2,el;

    el=end-start+1;


    if(el<2){
        return;
    }

    if(el==2){

        if(in.get_data(start)>in.get_data(end)){
            nn=in.get_data(start);
            in.set(start,in.get_data(end));
            in.set(end,nn);

            i1=dexes.get_data(start);
            dexes.set(start,dexes.get_data(end));
            dexes.set(end,i1);
        }

        return;
    }

    int i_mid,*i_use;

    i_mid=(start+end)/2;

    merge_sort(in,dexes,start,i_mid);
    merge_sort(in,dexes,i_mid+1,end);

    array_1d<T> buffer;
    array_1d<int> dex_buffer;

    buffer.set_where("merge_sort");
    dex_buffer.set_where("merge_sort");

    buffer.set_name("merge_sort_buffer");
    dex_buffer.set_name("merge_sort_dex_buffer");


    for(i1=start,i2=i_mid+1;i1<=i_mid || i2<=end;){

        if(i2>end){
            i_use=&i1;
        }
        else if(i1>i_mid){
            i_use=&i2;
        }
        else if(in.get_data(i1)<in.get_data(i2)){
            i_use=&i1;
        }
        else{
            i_use=&i2;
        }

        //printf("using %d -- %d %d -- %e %e\n",i_use[0],i1,i2,
        //in.get_data(i1),in.get_data(i2));

        buffer.add(in.get_data(i_use[0]));
        dex_buffer.add(dexes.get_data(i_use[0]));

        i_use[0]++;
    }

    //printf("start %d end %d buffer dim %d\n",start,end,buffer.get_dim());

    for(i1=0;i1<el;i1++){
        in.set(start+i1,buffer.get_data(i1));
        dexes.set(start+i1,dex_buffer.get_data(i1));
    }


}

template <typename T>
void sort(const array_1d<T> &in, array_1d<T> &sorted, array_1d<int> &dexes){
    sorted.set_dim(in.get_dim());
    int i;
    for(i=0;i<in.get_dim();i++){
        sorted.set(i,in.get_data(i));
    }
    merge_sort(sorted,dexes,0,in.get_dim()-1);
}

template <typename T>
double sort_and_check(const array_1d<T> &in, array_1d<T> &sorted, array_1d<int> &dexes){

    if(in.get_dim()!=dexes.get_dim()){
        printf("WARNING in sort_and_check in.dim %d dexes.dim %d\n",
        in.get_dim(),dexes.get_dim());

        exit(1);
    }

    if(in.get_dim()==0)return 0.0;

    in.set_where("sort_and_check");
    dexes.set_where("sort_and_check");
    in.set_where("sort_and_check");

    array_1d<int> dex_buffer;
    int i,j;

    dex_buffer.set_where("sort_and_check");

    dex_buffer.set_name("sort_and_check_dex_buffer");

    sorted.set_dim(in.get_dim());
    dex_buffer.set_dim(in.get_dim());

    for(i=0;i<in.get_dim();i++){
        dex_buffer.set(i,dexes.get_data(i));
        sorted.set(i,in.get_data(i));
    }

    merge_sort(sorted,dexes,0,in.get_dim()-1);

    double err,maxerr,aa,bb;

    int ifailure;

    for(i=0;i<in.get_dim();i++){
        if(i<in.get_dim()-1){
            if(sorted.get_data(i+1)<sorted.get_data(i)){
                printf("WARNING sort failed to get elements in proper order\n");

                ifailure=-1;

                throw ifailure;

            }
        }


        for(j=0;j<in.get_dim() && dexes.get_data(i)!=dex_buffer.get_data(j);j++);

        if(j==in.get_dim()){
            printf("WARNING could not find dex %d\n",dexes.get_data(i));
            ifailure=-1;
            throw ifailure;
        }

        if(dexes.get_data(i)!=dex_buffer.get_data(j)){
            printf("WARNING dexes did not line up %d %d\n",
            dexes.get_data(i),dex_buffer.get_data(j));

            ifailure=-1;

            throw ifailure;
        }

        aa=double(sorted.get_data(i));
        bb=double(in.get_data(j));

        err=fabs(aa-bb);
        if(fabs(aa)>0.0)err=err/fabs(aa);

        if(i==0 || err>maxerr){
            maxerr=err;
        }

    }

    if(maxerr>1.0e-12){
        printf("WARNING associative error in merge_sort was %e\n",maxerr);

        try{
          in.die(0);
        }
        catch(int iex){
           try{
               sorted.die(0);
           }
           catch(int jex){
               try{
                   dexes.die(0);
               }
               catch(int kex){

               };
           }
        }


        ifailure=-1;

        throw ifailure;
    }


    return maxerr;

}

template <typename T>
int get_dex(const array_1d<T> &xx, T target){
    int i;

    for(i=0;i<xx.get_dim()-1 && xx.get_data(i)<target;i++);

    if(i==0){
        if(target>xx.get_data(1)){
            printf("WARNING spuriously got zero in get_dex\n");
            printf("%e %e\n",double(target),double(xx.get_data(1)));
            exit(1);
        }
    }
    else if(i==xx.get_dim()-1){
        if(target<xx.get_data(i-1)){
            printf("WARNING spuriously got max in get_dex\n");
            printf("%e %e\n",double(target),double(xx.get_data(i-1)));
            exit(1);
        }
    }
    else if(target>xx.get_data(i) || target<xx.get_data(i-1)){
        printf("WARNING getdex failed %e -- %e %e\n",
        double(target),double(xx.get_data(i-1)),double(xx.get_data(i)));
        printf("%d %d\n",i,xx.get_dim()-1);
        exit(1);
    }


    if(i>0 && target<xx.get_data(i) && target-xx.get_data(i-1)<xx.get_data(i)-target){
        //printf("decrementing\n");
        i--;
    }

    return i;
}

template void merge_sort<double>(array_1d<double>&,array_1d<int>&,int,int);
template void merge_sort<int>(array_1d<int>&,array_1d<int>&,int,int);

template double sort_and_check<double>(const array_1d<double>&,array_1d<double>&,array_1d<int>&);
template double sort_and_check<int>(const array_1d<int>&,array_1d<int>&,array_1d<int>&);

template void sort<double>(const array_1d<double>&,array_1d<double>&,array_1d<int>&);
template void sort<int>(const array_1d<int>&,array_1d<int>&,array_1d<int>&);


template int get_dex(const array_1d<int>&,int);
template int get_dex(const array_1d<double>&,double);
