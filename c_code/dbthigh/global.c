int nseries;
int *series;
int newseries;
int newimage;
long ldum;  // seed for rand1()
int nbeans;
int firstb0done,iwantit;
int oops;
double threshlow,threshhigh,threshmask;
int dowater, dofat, dob0, dophase, dofwfrac, dowatfrac,rawimages, dicomimages, clobber, sloppy, reallysloppy, byimagenumber, fix, firstseries;
int verbose;
int debug;
int level;
char *thispatient;
char *basename;
double t1fat,t1wat,fwdiff;
double poff;
double pixelx,pixely;
char cmdline[1024];
char resultsdir[1024];

struct image {
  int series;
  int imageno;
  int rows;
  int cols;
  long accession;
  char *filename;
  char *sliceloc;
  char *imagepos;
  char *pixelspacing;
  char *slicethickness;
  char *acqtime;
  char *mediasopUID;
  char *sopUID;
  char *seriesUID;
  char *name;
  char *date;
  char imtype;
  float te;
  float tr;
  float flipangle;
  double larmor_freq;
  int echono;
  short *image;
  double *real;
  double *imag;
  double *b0;
  double *water;
  double *fat;
  int done;
  struct image *prev,*next;
};

struct tag {
  short group;
  short element;
  char VR[2];
};
struct image *basep;

char *seriesID[1024];  // why not, this is cheap

struct bean {
  double sliceloc;
  double b0mean;
  int order;
  int nvalidpts;
  double meanfwfrac;
  int milagro;  // this is the famous milagro bean field.
  int rows;
  int cols;
  struct bean *next;
  struct bean *prev;
};
struct bean *beanstar;

char destdir[1024];
