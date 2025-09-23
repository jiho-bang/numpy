#include <omp.h>
#include <x86intrin.h>

#include "compute.h"

// Computes the convolution of two matrices
int convolve(matrix_t *a_matrix, matrix_t *b_matrix, matrix_t **output_matrix) {
  // TODO: convolve matrix a and matrix b, and store the resulting matrix in
  // output_matrix
if (!output_matrix || !a_matrix || !b_matrix) {
        return -1;
    }
    int ac = a_matrix->cols;
    int ar = a_matrix->rows;
    int bc = b_matrix->cols;
    int br = b_matrix->rows;
    int* bd = b_matrix->data;
    int* ad = a_matrix->data;

  //if (!ad || !bd) {
      //return -1;
  //}

    matrix_t* b_flipped = malloc(sizeof(matrix_t));
    b_flipped->rows = bc;
    b_flipped->cols = br;
    b_flipped->data = malloc(sizeof(int) * br * bc);

    #pragma omp parallel for
    for (register int i = 0; i < br * bc; i++) {
      b_flipped->data[i] = bd[bc * br - i - 1];
    }

    int* bfd = b_flipped->data;

    // for (int i = 0; i < bc * br >> 1; i++) {
    // // for (int i = 0; i < (bc * br >> 1) >> 3 * 8; i+=8) {
    //   __m256i a = _mm256_loadu_si256((__m256i *) (b_matrix->data + i));
    //   __m256i b = _mm256_loadu_si256((__m256i *) (b_matrix->data + bc * br - i - 1));
    //   __m256i result = _mm256_permutevar8x32_epi32 (a, shuffle_pattern);
    //   __m256i result2 = _mm256_permutevar8x32_epi32 (b, shuffle_pattern);
    //   _mm256_storeu_si256((__m256i *)b_matrix->data + bc * br - i - 1, result);
    //   _mm256_storeu_si256((__m256i *)b_matrix->data + i, result2);
    // }
    // for (int i = (bc * br >> 1) >> 3 * 8; i < bc * br >> 1; i++) {
    //   int temp = bd[i];
    //   bd[i] = bd[bc * br - i - 1];
    //   bd[bc * br - i - 1] = temp;

    // //   __m256i a = _mm256_loadu_si256((__m256i *) (b_matrix->data + i));
    // //   __m256i b = _mm256_loadu_si256((__m256i *) (b_matrix->data + bc * br - i - 1));
    // //   __m256i result = _mm256_permutevar8x32_epi32 (a, shuffle_pattern);
    // //   __m256i result2 = _mm256_permutevar8x32_epi32 (b, shuffle_pattern);
    // //   _mm256_storeu_si256((__m256i *)b_matrix->data + bc * br - i - 1, result);
    // //   _mm256_storeu_si256((__m256i *)b_matrix->data + i, result2);
    // // }
    // // for (int i = (bc * br >> 1) >> 3 * 8; i < bc * br >> 1; i++) {
    // //   int temp = bd[i];
    // //   bd[i] = bd[bc * br - i - 1];
    // //   bd[bc * br - i - 1] = temp;
    // }


    int oc = ac - bc + 1; // output's col #
    int or = ar - br + 1;

    matrix_t* output = malloc(sizeof(matrix_t));
    output->rows = or;
    output->cols = oc;
    output->data = malloc(sizeof(int) * oc * or);


    #pragma omp parallel for           
    for (register int i = 0; i < or; i++) {
        #pragma omp parallel for 
        for (register int j = 0; j < oc; j++) {
           int result = 0;
            for (register int row = 0; row < br; row++) {
                __m256i sum_v = _mm256_setzero_si256();
                int col = 0;
                for (; col < ((bc >> 3)<< 3); col += 8) { 
                      __m256i a_data = _mm256_loadu_si256((__m256i *) (ad + (row+i) * ac + col+j));
                      __m256i b_data = _mm256_loadu_si256((__m256i *) (bfd + col + row * bc));
                      sum_v = _mm256_add_epi32(sum_v, _mm256_mullo_epi32(a_data, b_data));
                  }
                  int tmp_arr[8];
                  _mm256_storeu_si256((__m256i *)tmp_arr, sum_v);
                  tmp_arr[0] = (tmp_arr[0] + tmp_arr[1] + tmp_arr[2] + tmp_arr[3] + tmp_arr[4] + tmp_arr[5] + tmp_arr[6] + tmp_arr[7]);

                  // tailcase 
                  for (; col < bc; col++) {
                    tmp_arr[0] += ad[(row + i) * ac + j + col] * bfd[row * bc + col];
                  }
                  result += tmp_arr[0];
            }
            output->data[i * oc + j] = result;
        }
    }
    free(b_flipped->data);
    free(b_flipped);
    *output_matrix = output;  
    return 0;
}

// Executes a task
int execute_task(task_t *task) {
  matrix_t *a_matrix, *b_matrix, *output_matrix;

  char *a_matrix_path = get_a_matrix_path(task);
  if (read_matrix(a_matrix_path, &a_matrix)) {
    printf("Error reading matrix from %s\n", a_matrix_path);
    return -1;
  }
  free(a_matrix_path);

  char *b_matrix_path = get_b_matrix_path(task);
  if (read_matrix(b_matrix_path, &b_matrix)) {
    printf("Error reading matrix from %s\n", b_matrix_path);
    return -1;
  }
  free(b_matrix_path);

  if (convolve(a_matrix, b_matrix, &output_matrix)) {
    printf("convolve returned a non-zero integer\n");
    return -1;
  }

  char *output_matrix_path = get_output_matrix_path(task);
  if (write_matrix(output_matrix_path, output_matrix)) {
    printf("Error writing matrix to %s\n", output_matrix_path);
    return -1;
  }
  free(output_matrix_path);

  free(a_matrix->data);
  free(b_matrix->data);
  free(output_matrix->data);
  free(a_matrix);
  free(b_matrix);
  free(output_matrix);
  return 0;
}
