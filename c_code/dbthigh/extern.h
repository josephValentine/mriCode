#define FWDIFF 3.25 // default difference in resonant frequency of fat and water in ppm.

extern int nseries;
extern int *series;
extern int newseries;
extern int newimage;
extern long ldum;
extern int nbeans;
extern int firstb0done,iwantit;
extern int oops;
extern double threshlow,threshhigh,threshmask;
extern int dowater, dofat, dob0, dophase, dofwfrac, dowatfrac, rawimages, dicomimages, clobber, sloppy, reallysloppy, byimagenumber, fix, firstseries;
extern int verbose;
extern int debug;
extern int level;
extern char *thispatient;
extern char *basename;
extern double t1fat,t1wat,fwdiff;
extern double poff;
extern double pixelx,pixely;
extern char cmdline[1024];
extern char resultsdir[1024];

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
  unsigned short group;
  unsigned short element;
  char VR[2];
};
extern struct image *basep;

extern char *seriesID[1024];

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
extern struct bean *beanstar;

extern char destdir[1024];
extern char patname[512];
extern char sdate[32];
