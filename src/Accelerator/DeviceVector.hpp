// #include "cuda_error.h"
#ifndef LSMS_DEVICE_VECTOR_HPP
#define LSMS_DEVICE_VECTOR_HPP

#include "DeviceInterfaceCudaHip.hpp"
#include "deviceCheckError.hpp"
#include <vector>

  template <class T>
  class DeviceVector {
    public:
      typedef size_t       size_type;
      __inline__ DeviceVector() : N(0), data(0), owner(0) {}

      __inline__ DeviceVector( std::vector<T>& in): data(0) {
        *this=in;
      }
      __inline__ DeviceVector( size_type size ) {
        allocate(size*sizeof(T));
      }

      __inline__ ~DeviceVector() {
        free(); 
      }

      __inline__ __device__ T& operator[](size_type i) {
        return data[i];
      }

      __inline__ DeviceVector<T> &operator=(const std::vector<T>& in) {
        copy_async(in,0);
        return *this;
      }

      __inline__ void copy(const std::vector<T>& in) {
        size_type N=in.size(); 
        size_type num_bytes=N*sizeof(T);
        if(this->N!=N) {
          this->N=N;
          allocate(num_bytes);
        }
        deviceMemcpy(data,&in[0],num_bytes,deviceMemcpyHostToDevice);
        deviceCheckError();
      }
      
      __inline__ void copy_async(const std::vector<T>& in,deviceStream_t s) {
        size_type N=in.size(); 
        size_type num_bytes=N*sizeof(T);
        if(this->N!=N) {
          this->N=N;
          allocate(num_bytes);
        }
        deviceMemcpyAsync(data,&in[0],num_bytes,deviceMemcpyHostToDevice,s);
        deviceCheckError();
      }

      // We provide a few functions to return information about the vector
      __inline__ __host__ __device__ size_type size() const { return N; }
      __inline__ __host__ __device__ T* raw() const { return data; }

    private:
      __inline__ void allocate(size_type num_bytes) {
        if(num_bytes>0) {
          free();
          owner=this;
          deviceMalloc(&data,num_bytes);
          deviceCheckError();
        }
        else {     // DANGEROUS!! Might cause memory leak!
          owner=0;
          data=0;
          N=0;
        }
      }
      __inline__ void free() {
        if(owner==this && data!=0) {
          deviceFree(data);
          deviceCheckError();
          owner=0;
          data=0;
          N=0;
        }
      }

      DeviceVector<T> *owner;
      size_type N;
      T* data;
  };

#endif

