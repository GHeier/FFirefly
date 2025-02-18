#ifndef CONDENSEDMATTERFIELD_H
#define CONDENSEDMATTERFIELD_H

#ifdef __cplusplus
extern "C" {
#endif

// Initialize Julia runtime (important when calling from C++)
void julia_init();

// Create CMF object from file
void *create_CMF_from_file(const char *filename);

// Evaluate CMF at given (q, w)
double eval_CMF(void *cmf_obj, double *q, int q_length, double w);

// Free CMF object memory
void free_CMF(void *cmf_obj);

#ifdef __cplusplus
}
#endif

#endif /* CONDENSEDMATTERFIELD_H */

