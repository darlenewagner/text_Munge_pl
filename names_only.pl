use strict;

## Removes numerical digits and non-grammatical special characters

while(<STDIN>)
{
    $_ =~ s/(\d|\†|\‡|\§|\*)//g;
    $_ =~ s/,,/,/g;
    print;

}
