/* Make one declaration of the init function for each model. */
void NLP_model_init(struct s_model *model, FILE *fptr);
void NLP_u_model_init(struct s_model *model, FILE *fptr);
void PHARM_model_init(struct s_model *model, FILE *fptr);

/* Add one line to all_models giving the model name and address to the
 * initialization function. */
struct s_all_models {
    char *name;
    void (*init_fcn)(struct s_model *model, FILE *fptr);
} all_models[] = {
    {"NLP", &NLP_model_init },
    {"NLP_u", &NLP_u_model_init },
    {"PHARM", &PHARM_model_init }
};

#define NUM_MODELS (sizeof all_models / sizeof(struct s_all_models))
