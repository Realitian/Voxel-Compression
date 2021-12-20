#ifndef SL_GETOPT_HPP
#define SL_GETOPT_HPP

namespace sl {

  extern int opterr;		/* if error message should be printed */
  extern int optind;		/* index into parent argv vector */
  extern int optopt;		/* character checked for validity */
  extern int optreset;		/* reset getopt */
  extern char *optarg;		/* argument associated with option */
  
  struct option {
    const char *name;
    int has_arg;
    int *flag;
    int val;
  };
  
#define sl_no_argument       0
#define sl_required_argument 1
#define sl_optional_argument 2

  extern int getopt(int argc, char * const argv[],
		    const char *optstring);
  extern int getopt_long(int argc, char * const argv[],
			 const char *optstring,
			 const struct option *longopts, int *longindex);

} // namespace sl 

#endif /* SL_GETOPT_HPP */
