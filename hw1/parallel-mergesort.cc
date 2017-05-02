#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>

#include "sort.hh"

void mergeSort(keytype* A, keytype* B, int start, int end, int thresold);
int comp(const void* a, const void* b);
void pMerge(keytype* A, keytype* B, int start1, int end1, int start2, int end2, int start, int thresold);
void sMerge(keytype* A, keytype* B, int start1, int end1, int start2, int end2, int start);

void parallelSort (int N, keytype* A) {
    keytype* B = new keytype[N];
    
    #pragma omp parallel
    {
        int noThreads = omp_get_num_threads();
        
    #pragma omp single
        {
            mergeSort(A, B, 0, N-1, N/noThreads);
        }
    }
}

void mergeSort(keytype* A, keytype* B, int start, int end, int thresold) {
    
    int n = end - start + 1;
    
    if ( n <= thresold) {
        qsort(A + start, n, sizeof(keytype), comp);
        return;
    } else {
        int mid = start + (end - start) / 2;
        
        #pragma omp task
        {
            mergeSort(A, B, start, mid, thresold);
        }
        mergeSort(A, B, mid + 1, end, thresold);
           
        #pragma omp taskwait
        pMerge(A, B, start, mid, mid + 1, end, start, thresold);
        
        memcpy(A + start, B + start, (end - start + 1) * sizeof(keytype));
    }
    
}

int comp(const void* a, const void* b)
{
    return ( *(int*)a - *(int*)b );
}

int bSearch(keytype key, keytype* arr, int start, int end) {
    
    int low = start;
    int high = start > (end+1) ? start : (end+1);
    int mid;
    
    while (low < high) {
        
        mid = low + (high - low) / 2;
        
        if (key <= arr[mid]) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    return high;
}

void pMerge(keytype* A,  keytype* B, int start1, int end1, int start2, int end2, int start, int thresold) {
    
    int n1 = end1 - start1 + 1;
    int n2 = end2 - start2 + 1;
    
    if((n1+n2) > thresold){
        
        if (n1 < n2) {
            
            int temp;
            
            temp = start1;
            start1 = start2;
            start2 = temp;
            
            temp = end1;
            end1 = end2;
            end2 = temp;
            
            temp = n1;
            n1 = n2;
            n2 = temp;
            
        }
        
        if (n1 == 0) {
            
            return;
            
        } else {
            
            int mid = (start1 + end1) / 2;
            int idx = bSearch(A[mid], A, start2, end2);
            int q = start + (mid - start1) + (idx - start2);
            B[q] = A[mid];
            
            #pragma omp task
            {
                pMerge(A, B, start1, mid - 1, start2, idx - 1, start, thresold);
            }
            pMerge(A, B, mid + 1, end1, idx, end2, q + 1, thresold); 
            
            #pragma omp taskwait
        }
    } else {
        sMerge(A, B, start1, end1, start2, end2, start);
    } 
    
}


void sMerge(keytype* A, keytype* B, int start1, int end1, int start2, int end2, int start) {
    
    int i = start1, j = start2, k = start;
    
    while ( (i <= end1) && (j <= end2) ) {
        
        if ( A[i] <= A[j] ) {
            B[k++] = A[i++];
        } else {
            B[k++] = A[j++];
        }
    }
    
    while (i <= end1) {
        B[k++] = A[i++];
    }
    while (j <= end2) {
        B[k++] = A[j++];
    }
}
