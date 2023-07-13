#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


double left[500000];
double right[500000];

void merge(double *a, int s, int m, int e)
{
    int n1 = m - s + 1;
    int n2 = e - m;

    

    for (int i = 0; i < n1; i++)
        left[i] = a[s + i];
    for (int i = 0; i < n2; i++)
        right[i] = a[i + m + 1];

    int i1 = 0, i2 = 0, i = s;
    while (i1 < n1 && i2 < n2)
    {
        if (left[i1] <= right[i2])
        {
            a[i++] = left[i1++];
        }
        else
        {
            a[i++] = right[i2++];
        }
    }
    while (i1 < n1)
        a[i++] = left[i1++];
    while (i2 < n2)
        a[i++] = right[i2++];
}

double minima(double a, double b)
{
    if (a < b)
        return a;
    else
        return b;
}
double maxima(double a, double b)
{
    if (a > b)
        return a;
    else
        return b;
}

void insertionSort(double *arr, int start,int end)
{
    int i, j;
    double key;
    for (i = start; i <= end; i++) {
        key = arr[i];
        j = i - 1;
 
        while (j >= start && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void ImprovedMergesort(double *ar, int start, int end, int n)
{

    if(start+3<end)//Use insertion sort for smaller sized array.
    {
        insertionSort(ar,start,end);
    }
   if (start < end)
        {
            int mid = (start + end) / 2;
            ImprovedMergesort(ar, start, mid, n);
            ImprovedMergesort(ar, mid + 1, end, n);
            if (ar[mid + 1] < ar[mid])
            merge(ar, start, mid, end);
        }
    }


    int main()
    {
        srand((unsigned)time(NULL));
        int MAX_RANGE = 10000;
        FILE *fptr;
        fptr = fopen("ImprovedMergeSortResults.txt", "w+");
        int n;
        for (n = 100; n <= 100; n *= 10)
        {

            fprintf(fptr, "\n Results for value of n = %d \n", n);
            fprintf(fptr, " 2n loge n : %f \n", 2 * n * log10(n) / log10(exp(1)));
            fprintf(fptr, " n log2 n : %f \n", n * log2(n));

            if (fptr == NULL)
            {
                printf("Error !File doesn't exist");
                exit(1);
            }

            clock_t t;
            clock_t time_taken;
            t = clock();
            long long int total_comparisons = 0;

            for (int k = 1; k <= 1; k++)
            {
                double *arr = (double *)malloc(sizeof(double) * n);

                for (int i = 0; i < n; i++)
                {
                    double random_num = ((double)rand() / (double)RAND_MAX) * MAX_RANGE;
                    arr[i] = random_num;
                }
                // for(int i=0;i<n;i++)
                // printf("%f ",arr[i]);

                ImprovedMergesort(arr, 0, n - 1, n);
                for(int i=0;i<n;i++)
                printf("%f ",arr[i]);
                free(arr);

                
            }

            time_taken = (clock() - t);
            int time_taken_microseconds = ((double)time_taken * 1000000) / (double)CLOCKS_PER_SEC;

            fprintf(fptr, "\n\nTime Taken : %f μs\n", time_taken_microseconds / 500.0);
        }
    }
