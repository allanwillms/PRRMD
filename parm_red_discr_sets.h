/* Make one declaration of the init function for each discretization set. */
void A1OUT_discr_set_init(struct s_discr_set *discr_set, FILE *fptr);

/* Add one line to all_discr_sets giving the discretization set name and address 
 * to the initialization function. */
struct s_all_discr_sets {
    char *name;
    void (*init_fcn)(struct s_discr_set *discr_set, FILE *fptr);
} all_discr_sets[] = {
    {"A1OUT", &A1OUT_discr_set_init }
};

#define NUM_DSETS (sizeof all_discr_sets / sizeof(struct s_all_discr_sets))
