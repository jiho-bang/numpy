#include "compute.h"

// Computes the convolution of two matrices
int convolve(matrix_t *a_matrix, matrix_t *b_matrix, matrix_t **output_matrix) {
  // TODO: convolve matrix a and matrix b, and store the resulting matrix in
  // output_matrix
  if (!output_matrix || !a_matrix || !b_matrix) {
      return -1;
  }
  
  uint32_t ac = a_matrix->cols;
  uint32_t ar = a_matrix->rows;
  uint32_t bc = b_matrix->cols;
  uint32_t br = b_matrix->rows;

  if (!ac || !ar || !bc || !br) {
      return -1;
  }

  int32_t* ad = a_matrix->data;
  int32_t* bd = b_matrix->data;

  if (!ad || !bd) {
      return -1;
  }

  for (uint32_t i = 0; i < br; i++) {
      for (uint32_t j = 0; j < bc >> 1; j++) {
          int32_t temp = bd[i * bc + j];
          bd[i * bc + j] = bd[(i + 1) * bc - j - 1];
          bd[(i + 1) * bc - j - 1] = temp;
       }
  }
  for (uint32_t i = 0; i < bc; i++) {
      for (uint32_t j = 0; j < br >> 1; j++) {
          int32_t temp = bd[j * bc + i];
          bd[j * bc + i] = bd[(br - j - 1) * bc + i];
          bd[(br - j - 1) * bc + i] = temp;
      }
  }

  matrix_t* output = malloc(sizeof(matrix_t));
  output->rows = ar - br + 1;
  output->cols = ac - bc + 1;
  output->data = malloc(sizeof(int32_t) * output->rows * output->cols);
  int32_t* temp = malloc(sizeof(int32_t) * br * bc);
  int n = br * bc;

  for (uint32_t i = 0, x = 0; i < output->rows; i++) {
      for (uint32_t j = 0; j < output->cols; j++) {
          for (uint32_t k = i, y = 0; k < i + br; k++) {
              for (uint32_t l = j; l < j + bc; l++) {
                  temp[y++] = ad[k * ac + l];
              }
          }
          int32_t dp = 0;
          for (uint32_t d = 0; d < n; d++) {
              dp += temp[d] * bd[d];
          }
          output->data[x++] = dp;
      }
  }
  free(temp);
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
