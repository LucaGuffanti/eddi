#include "eddi_density_readers.h"
#include "eddi_density_writers.h"

int main(int argc, char** argv)
{
    eddi_density_field_t density_field;
    eddi_molecule_t molecule;

    char* binary_name = malloc(strlen(argv[1]) + 5);

    int i = 0;
    for (; argv[1][i+5]!='\0'; ++i)
    {
        binary_name[i] = argv[1][i];
    }

    binary_name[i++] = '.';
    binary_name[i++] = 'b';
    binary_name[i++] = 'i';
    binary_name[i++] = 'n';
    binary_name[i] = '\0';

    printf("%s -> %s\n", argv[1], binary_name);
    eddi_read_gaussian_cube(argv[1], &density_field, &molecule);
    eddi_write_binary(binary_name, &density_field);

    return 0;
}