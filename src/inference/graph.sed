#!/usr/bin/sed -f

# This program scans a C++ file and looks for functions b(...) that are called
# inside the body of another function a(...) { ... b(...) }. Each instance of b
# is replaced with 'a -> b': so that
#    returntype parent(...) { ... child(...) } --> 'parent -> child'

# match function definitions:
# the pattern to find is '^type name(' and we want to get the name:
# match the line that defines the function:
/^[a-z].* [a-z_-]\+(/{
    # replace '^type name(' with 'name':
    s/^[a-z].* \([a-z_-]\+\)(.*/\1/;
    h; # put the name of the parent in hold space
    :loop; # start a loop...
    # If a child function is called on this line:
    /[a-z](/{
        # find the name of the child and delete everything else:
        s/.*[ :-]\([a-z_-]\+\)(.*$/\1/;
        # glue the parent name onto the end of the child:
        G;
        # remove a newline:
        s/\n/ /;
        # rewrite the string to 'parent -> child':
        s/\(.*\) \(.*\)/\2 -> \1/;
        # print the result:
        p;
    };
    # read in the next line:
    n;
# when we find a close paren on its own, we have reached the end of the function
# and need to exit the loop. so, as long as we have NOT found the exit pattern,
# repeat the loop.
    /^}$/!b loop; # loop while we are in the function.
}
