#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

long long int comparisons; 

//These two arrays have been created for Imoproved Merge Sort.
double leftI[500000];  // Two global arrays to copy left and right
double rightI[500000]; // sorted parts of the main array 

void mergeImproved(double *a, int s, int m, int e) // Merge Function for improved merge sort
{
    int n1 = m - s + 1;
    int n2 = e - m;

    for (int i = 0; i < n1; i++)
        leftI[i] = a[s + i];
    for (int i = 0; i < n2; i++)
        rightI[i] = a[i + m + 1];

    int i1 = 0, i2 = 0, i = s;
    while (i1 < n1 && i2 < n2)
    {
        if (leftI[i1] <= rightI[i2])
        {
            a[i++] = leftI[i1++];
        }
        else
        {
            a[i++] = rightI[i2++];
        }
    }
    while(i1<n1)
    a[i++] = leftI[i1++];
    while(i2<n2)
    a[i++] = rightI[i2++];
}

void merge(double *a, int s, int m, int e) // Merge Function for merge sort
{
    int n1 = m - s + 1;
    int n2 = e - m;

    double  left[n1]; // Creating array multiple times for trivial merge sort
    double right[n2];

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
        comparisons++; // Keep increasing the number of comparisons after each iteration
    }
    while (i1 < n1)
        a[i++] = left[i1++];
    while (i2 < n2)
        a[i++] = right[i2++];
}

void swap(double *x, double *y) // Helper function for swapping two variables
{
    double temp = *x;
    *x = *y;
    *y = temp;
}

int partition(double *arr, int l, int r) // Optimized function for partition around a pivot (the first element of the array)
{                                        // element of the array) in this case
    double pivot = arr[l];
    int i = l + 1;
    int end = r;

    while (i <= end)
    {
        if (arr[i] > pivot)
            swap(&arr[end--], &arr[i]);
        else
            i++;

        comparisons++;
    }

    swap(&arr[l], &arr[end]);
    return end;
}

// ***MERGE SORT***
void MergeSort(double *arr, int l, int r)
{
    if (l < r)
    {
        int mid = (l + r) / 2;
        MergeSort(arr, l, mid);
        MergeSort(arr, mid + 1, r);
        
        merge(arr, l, mid, r);
    }
}

// ***QUICK SORT***
void QuickSort(double *arr, int l, int r)
{

    if (l < r)
    {
        int i = partition(arr, l, r);
        
        QuickSort(arr, l, i - 1);
        QuickSort(arr, i + 1, r);
    }
}


// ***IMPROVED MERGE SORT***

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

void ImprovedMergeSort(double* arr, int l, int r)
{
    if(l+3<r)//Use insertion sort for smaller sized array.
    {
        insertionSort(arr,l,r);
    }
    if(l<r)
    {
        int mid = (l+r)/2;
        ImprovedMergeSort(arr, l, mid);
        ImprovedMergeSort(arr, mid+1, r);
        
        if(arr[mid]>arr[mid+1])//This step reduces the redudant memory overheads.
        mergeImproved(arr, l, mid,r); 
    }
}
int main()
{

    srand((unsigned)time(NULL));
    int MAX_RANGE = 10000;
    FILE *fptr;
    fptr = fopen("AssignmentSortResults.txt", "w+"); // To print the output in a text file
    int n;
    for (n = 100; n <= 1000000; n *= 10)
    {
        fprintf(fptr, "\n Results for value of n = %d \n", n);

        if (fptr == NULL)
        {
            printf("Error !File doesn't exist");
            exit(1);
        }
        int pos=0;
        // Creating arrays to store time taken by the 3 algorithms for each iteration
        double *checkerquick = (double *)malloc(sizeof(double) * 500); 
        double *checkermerge = (double *)malloc(sizeof(double) * 500);
        double *checkerimpmerge = (double *)malloc(sizeof(double) * 500);
        clock_t t;
        clock_t time_taken;
        long long int total_comparisons_quick = 0;
        long long int total_comparisons_merge = 0;

        for (int k = 1; k <= 500; k++)
        {

            double *arrquick = (double *)malloc(sizeof(double) * n);
            double *arrmerge = (double *)malloc(sizeof(double) * n);
            double *arrimpmerge = (double *)malloc(sizeof(double) * n);

            //Creating all the arrays.
            for (int i = 0; i < n; i++)
            {
                double random_num = ((double)rand() / (double)RAND_MAX) * MAX_RANGE;
                arrquick[i] = random_num;
                arrmerge[i] = random_num;
                arrimpmerge[i]=random_num;
            }

            // Doing all things for quick sort.
            
            t = clock();
            comparisons = 0;
            QuickSort(arrquick, 0, n - 1);
            time_taken = (clock() - t);
            total_comparisons_quick += comparisons;
            int time_taken_microseconds_quick = ((double)time_taken * 1000000) / (double)CLOCKS_PER_SEC;
            checkerquick[pos] = time_taken_microseconds_quick;

            // Now we do all the things for merge sort.
            
            t = clock();
            comparisons = 0;
            MergeSort(arrmerge, 0, n - 1);
            time_taken = (clock() - t);
            total_comparisons_merge += comparisons;
            int time_taken_microseconds_merge = ((double)time_taken * 1000000) / (double)CLOCKS_PER_SEC;
            checkermerge[pos] = time_taken_microseconds_merge;

            // Now we do all the things for Improved merge sort.
            
            t = clock();
            ImprovedMergeSort(arrimpmerge, 0, n - 1);
            time_taken = (clock() - t);
            int time_taken_microseconds_Impmerge = ((double)time_taken * 1000000) / (double)CLOCKS_PER_SEC;
            checkerimpmerge[pos] = time_taken_microseconds_Impmerge;

            //free all the arrays that were created.
            pos++;
            free(arrquick);
            free(arrmerge);
            free(arrimpmerge);
        }

        fprintf(fptr, " 2n loge n : %f \n", 2 * n * log10(n) / log10(exp(1)));
        fprintf(fptr, " n log2 n : %f \n", n * log2(n));

        // Finding average running time for all the sorting techniques.
        double totaltimequick = 0;
        for (int j = 0; j < 500; j++)
            totaltimequick += checkerquick[j];
        double avgtimequick = totaltimequick / 500.0;

        double totaltimemerge = 0;
        for (int j = 0; j < 500; j++)
            totaltimemerge += checkermerge[j];
        double avgtimemerge = totaltimemerge / 500.0;

        double totaltimeimpmerge = 0;
        for (int j = 0; j < 500; j++)
            totaltimeimpmerge += checkerimpmerge[j];
        double avgtimeimpmerge = totaltimeimpmerge / 500.0;

        fprintf(fptr, "Average number of comparisons made Quick: %lld\n\n ", total_comparisons_quick / 500);
        fprintf(fptr, "Average number of comparisons made Merge: %lld\n\n ", total_comparisons_merge / 500);

        fprintf(fptr, "\n\nAverage Running time for Quick Sort : %f μs\n", avgtimequick);
        fprintf(fptr, "\n\nAverage Running time for Merge Sort : %f μs\n", avgtimemerge);
        fprintf(fptr, "\n\nAverage Running time for Improved Merge Sort : %f μs\n", avgtimeimpmerge);

        //Calculating the number of times mergesort outperformed quick sort.
        int mergewin=0;
        for(int i=0;i<500;i++)
        {
            if(checkermerge[i]<checkerquick[i])
            mergewin++;
        }

        //Calculating the number of times mergesort outperformed quick sort.
        int mergeimpwin=0;
        for(int i=0;i<500;i++)
        {
            if(checkerimpmerge[i]<checkerquick[i])
            mergeimpwin++;
        }

        fprintf(fptr,"\n\n Merge wins for %d times \n\n",mergewin);
        fprintf(fptr,"\n\n Improved Merge wins for %d times \n\n",mergeimpwin);

        int percent5 = 0, percent10 = 0, percent20 = 0, percent30 = 0, percent50 = 0, percent100 = 0;

        for (int j = 0; j < 500; j++)
        {
            if (checkerquick[j] > 2 * avgtimequick)
                percent100++;
            if (checkerquick[j] > 1.5 * avgtimequick)
                percent50++;
            if (checkerquick[j] > 1.3 * avgtimequick)
                percent30++;
            if (checkerquick[j] > 1.2 * avgtimequick)
                percent20++;
            if (checkerquick[j] > 1.1 * avgtimequick)
                percent10++;
            if (checkerquick[j] > 1.05 * avgtimequick)
                percent5++;
        }
        // free(checker);
        fprintf(fptr, "No. of cases where run time exceeds average by 5 %d\n", percent5);
        fprintf(fptr, "No. of cases where run time exceeds average by 10 %d\n", percent10);
        fprintf(fptr, "No. of cases where run time exceeds average by 20 %d\n", percent20);
        fprintf(fptr, "No. of cases where run time exceeds average by 30 %d\n", percent30);
        fprintf(fptr, "No. of cases where run time exceeds average by 50 %d\n", percent50);
        fprintf(fptr, "No. of cases where run time exceeds average by 100 %d\n", percent100);
    }

    printf("\n\n\n\n");

    printf("Answers for 1.2 part: \n");
    for (n = 100000; n <= 900000; n += 200000)
    {

        int avg = 0;
        double *arrquick = (double *)malloc(sizeof(double) * n);
        // Creating all the arrays.
        for (int k = 1; k <= 500; k++)
        {
            for (int i = 0; i < n; i++)
            {
                double random_num = ((double)rand() / (double)RAND_MAX) * MAX_RANGE;
                arrquick[i] = random_num;
            }

            // Doing all things for quick sort.
            clock_t t;
            clock_t time_taken;
            t = clock();
            comparisons = 0;
            QuickSort(arrquick, 0, n - 1);
            time_taken = (clock() - t);
            int time_taken_microseconds_quick = ((double)time_taken * 1000000) / (double)CLOCKS_PER_SEC;
            avg += time_taken_microseconds_quick;

            
        }
            printf("Average running time for %d : %f\n ",n,avg/500.0);
            printf("nlogn = %f \n\n",n*log10(n));
    }
}