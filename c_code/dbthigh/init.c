#include<stdlib.h>
#include<stdio.h>
#include<sys/types.h>
#include<dirent.h>
#include<string.h>
#include "extern.h"
#define WATER_T1 640.0  // T1 of water component of leg muscle in ms, measured from Smith's leg
// ("What's the name of his other leg?")
#define FAT_T1 340.0    // T1 of fat component of leg in ms.

double atof();

int dicomreadloopdone;
int iswantedloopdone;

init(argc,argv)
     int argc;
     char **argv;
{
  int i,j,k;
  DIR *dirp;
  struct dirent *direntp;
  FILE *fp;
  char fourchar[5];
  char filename[PATH_MAX + 1];
  struct image *ip;
  int gotseries;
  double atof();

  level = 0;
  cmdline[0] = '\0';
  sprintf(cmdline,"%s",argv[0]);
  for(i=1; i<argc; i++)
    sprintf(cmdline,"%s %s",cmdline,argv[i]);
  printf("command line: %s\n",cmdline);

  destdir[0] = '\0';
  threshlow = 5.0;   // threshold values as percent of max
  threshhigh = 20.0;
  threshmask = -1.0;
  ldum = -1;
  basep = (struct image *) 0;
  oops = 0;
  gotseries = 0;
  sloppy = reallysloppy = byimagenumber = fix = clobber = dowater = dofat = dob0 = dophase = dofwfrac = dowatfrac = rawimages = dicomimages = debug = 0;
  basename = (char *)0;
  t1fat = t1wat = poff = 0.0;
  pixelx = pixely = 0.0;
  verbose = 0;
  resultsdir[0] = '\0';
  firstseries = 100;
  fwdiff = FWDIFF;
  for(i=0; i<6; i++)
    seriesID[i] = (char *)0;

  if(argc < 3)
    usage(argv);
  for(i=1; i<argc; i++) {
    if(!strcmp(argv[i],"-v")) {
      verbose++;
      continue;
    }
    if(!strcmp(argv[i],"-poff")) {
      i++;
      poff = atof(argv[i]);
      continue;
    }
    if(!strcmp(argv[i],"-dir")) {
      i++;
      strcpy(destdir,argv[i]);
      continue;
    }
    if(!strcmp(argv[i],"-fwdiff")) {
      i++;
      fwdiff = atof(argv[i]);
      continue;
    }
    if(!strcmp(argv[i],"-oops")) {
      i++;
      oops = atoi(argv[i]);
      printf("oops = %d\n",oops);
      continue;
    }
    if(!strcmp(argv[i],"-sloppy")) {
      sloppy++;
      continue;
    }
    if(!strcmp(argv[i],"-debug")) {
      debug++;
      continue;
    }
    if(!strcmp(argv[i],"-reallysloppy")) {
      reallysloppy++;
      continue;
    }
    if(!strcmp(argv[i],"-byno")) {
      byimagenumber++;
      continue;
    }
    if(!strcmp(argv[i],"-fix")) {
      fix++;
      continue;
    }
    if(!strcmp(argv[i],"-mask")) {
      i++;
      threshmask = atof(argv[i]);
      continue;
    }
    if(!strcmp(argv[i],"-name")) {
      basename = malloc(strlen(argv[++i])+1);
      strcpy(basename,argv[i]);
      continue;
    }
    if(!strcmp(argv[i],"-fw")) {
      dofwfrac++;
      continue;
    }
    if(!strcmp(argv[i],"-wf")) {
      dowatfrac++;
      continue;
    }
    if(!strcmp(argv[i],"-w")) {
      dowater++;
      continue;
    }
    if(!strcmp(argv[i],"-f")) {
      dofat++;
      continue;
    }
    if(!strcmp(argv[i],"-b")) {
      dob0++;
      continue;
    }
    if(!strcmp(argv[i],"-p")) {
      dophase++;
      continue;
    }
    if(!strcmp(argv[i],"-d")) {
      dicomimages++;
      continue;
    }
    if(!strcmp(argv[i],"-s")) {
      rawimages++;
      continue;
    }
    if(!strcmp(argv[i],"-c")) {
      clobber++;
      printf("clobbering prior output files with reckless abandon...\n");
      continue;
    }
    if(!strcmp(argv[i],"-n")) {
      i++;
      firstseries = atoi(argv[i]);
      printf("first addition series number is %d",firstseries);
      if(firstseries < 0) {
	fprintf(stderr,"negative series numbers not allowed, bailing...\n");
	exit(1);
      }
      continue;
    }
    if(!strcmp(argv[i],"-t1fat")) {
      i++;
      t1fat = atof(argv[i]);
      printf("fat t1 specified as %g ms\n",t1fat);
      continue;
    }
    if(!strcmp(argv[i],"-t1wat")) {
      i++;
      t1wat = atof(argv[i]);
      printf("water t1 specified as %g ms\n",t1wat);
      continue;
    }
    if(!strcmp(argv[i],"-lt")) {
      i++;
      threshlow = atof(argv[i]);
      printf("low threshold set to %g\%\n",threshlow);
      continue;
    }
    if(!strcmp(argv[i],"-ht")) {
      i++;
      threshhigh = atof(argv[i]);
      printf("high threshold set to %g\%\n",threshhigh);
      continue;
    }
    if(!gotseries) {
      getseries(argv[i]);
      gotseries++;
      continue;
    }
    if((dirp = opendir(argv[i])) == NULL) {
      fprintf(stderr,"could not open directory %s: ",argv[i]);
      perror("");
      exit(1);
    }
    printf("directory %s:\n",argv[i]);
    for(direntp = readdir(dirp); direntp != NULL; direntp = readdir(dirp)) {
      // files on CDROM seem to be type "unknown"
      if((direntp->d_type != DT_REG) && (direntp->d_type != DT_UNKNOWN))
	continue;
      sprintf(filename,"%s/%s",argv[i],direntp->d_name);
      if((fp = fopen(filename,"r")) == NULL) {
	fprintf(stderr,"could not open file %s: ",filename);
	perror("");
	exit(1);
      }
      if(fseek(fp,128L,SEEK_CUR) == -1) {
	perror("fseek error");
	exit(1);
      }
      if(fread(fourchar,1,4,fp) != 4) {
	// this will fail if we have seeked beyond the end of file, i.e. files of length < 128 bytes
	continue;
      }
      fourchar[4] = '\0';
      if(strcmp(fourchar,"DICM")) {
	continue;
      }
      if(!iswanted(fp)) {
	fclose(fp);
	continue;
      }

      if(basep == (struct image *)0) {
	if((basep = (struct image *)malloc(sizeof(struct image))) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	basep->prev = (struct image *)0;
	basep->next = (struct image *)0;
	ip = basep;
      }
      else {
	if((ip->next = (struct image *)malloc(sizeof(struct image))) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	ip->next->prev = ip;
	ip = ip->next;
	ip->next = (struct image *)0;
      }

      ip->done = 0;
      readdicom(ip,fp);
      if((ip->filename = malloc(strlen(filename) + 1)) == NULL) {
	perror("malloc failed");
	exit(1);
      }
      strcpy(ip->filename,filename);

      fclose(fp);
    }
    closedir(dirp);
  }
  printf("patient name: %s\n",ip->name);
  if((thispatient = malloc(strlen(ip->name)+1)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  strcpy(thispatient,ip->name);
  thispatient[strlen(ip->name)] = '\0';
  if(!dofat && !dowater && !dob0 && !dowatfrac && !dofwfrac)
    dofwfrac++;
  if(!dicomimages && !rawimages)
    dicomimages++;
  if(t1wat == 0.0)
    t1wat = WATER_T1;
  if(t1fat == 0.0) 
    t1fat = FAT_T1;
}

usage(argv) 
     char **argv;
{
  fprintf(stderr,"usage: %s [-oops n] [-lt thresh] [-ht thresh] [-fw] [-f] [-w] [-b] [-p] [-d] [-s] [other options] [-dir destdir] [-name bname] series dir1 [dir2 ...]\n",argv[0]);
  fprintf(stderr,"\t-oops n\toffset phase of center field map by n*2pi (i.e. we centered on fat instead of water frequency\n");
  fprintf(stderr,"\t-lt\tlow threshold for phase unwrapping in %% of max [5.0]\n");
  fprintf(stderr,"\t-ht\thigh threshold for phase unwrapping in %% of max [20.0]\n");
  fprintf(stderr,"\t-fw\tcreate fat-water fraction images in units of 0.1\% (default)\n");
  fprintf(stderr,"\t-wf\tcreate water fraction images in units of 0.1\%\n");
  fprintf(stderr,"\t-f\tcreate fat-only images from Dixon method\n");
  fprintf(stderr,"\t-w\tcreate water-only images from Dixon method\n");
  fprintf(stderr,"\t-w\tcreate phase images for water or fat only images\n");
  fprintf(stderr,"\t-b\tcreate B0 images\n");
  fprintf(stderr,"\t-d\tcreate DICOM image files (default)\n");
  fprintf(stderr,"\t-s\tcreate raw data files for Scion image\n");
  fprintf(stderr,"\t-sloppy\treport center freq mismatches but keep going\n");
  fprintf(stderr,"\t-reallysloppy\tallow all sorts of things to not match\n");
  fprintf(stderr,"\t-byno\tmatch slices by image number instead of slice location\n");
  fprintf(stderr,"\t-fix\tadjust phase data for frequency mismatch\n");
  fprintf(stderr,"\t-dir destdir\tputs files into subdirectories of destdir by pat name and date\n");
  fprintf(stderr,"\t-c\tclobber results files silently\n");
  fprintf(stderr,"\t-debug\toutput tons of dicom tag information for debugging purposes\n");
  fprintf(stderr,"\t-n x\tmake first new series number x [100]\n");
  fprintf(stderr,"\t-name bname\tmake base name for output files bname []\n");
  fprintf(stderr,"\t-t1fat number\tuse number as t1 for fat resonance in ms\n");
  fprintf(stderr,"\t-t1wat number\t use number as t1 for water resonance in ms\n");
  fprintf(stderr,"\t-fwdiff number\t use number as freq diff between fat and water in ppm\n");
  fprintf(stderr,"\t-poff number\tadd poff degrees of phase to second echo image\n");
  fprintf(stderr,"\t-mask number\tset mag threshold in percent for volume fraction images\n");
  fprintf(stderr,"\t-v\tverbose\n");
  fprintf(stderr,"\tseries:\tcomma separated list of series numbers, e.g. 1,7,9 (no spaces)\n");
  fprintf(stderr,"\tdir1 [dir2...]:\tlist of directories to search\n");
  exit(1);
}

iswanted(fp)
     FILE *fp;
{
  fseek(fp,132L,SEEK_SET);
  iwantit = 0;
  iswantedloopdone = 0;
  iswantedloop(&fp,-1);
  return(iwantit);
}

iswantedloop(fpp,itemlength)
     FILE **fpp;
     long itemlength;
{
  long start;
  long sqstart,sqlength;
  struct tag {
    unsigned short group;
    unsigned short element;
  } mytag;
  char VR[2];
  unsigned int length;
  unsigned short slength;
  unsigned short sdum;
  unsigned int intdum;
  char cbuf[256];
  int i,item,thisseries;

  if(itemlength == 0)
    return;
  start = ftell(*fpp);
  for(; !iswantedloopdone;) {
    if(fread(&mytag,sizeof(struct tag),1,*fpp) != 1) {
	if(feof(*fpp))
	  return;
	perror("read error");
	exit(1);
      }

      if(mytag.group == 0xFFFE)
	VR[0] = VR[1] = '\0';
      else
	if(fread(&VR[0],sizeof(VR),1,*fpp) != 1) {
	  if(feof(*fpp))
	    return;
	  perror("read error");
	  exit(1);
	}

      if(!strncmp(VR,"OB",2) || !strncmp(VR,"OW",2) || !strncmp(VR,"OF",2)
	 || !strncmp(VR,"SQ",2) || !strncmp(VR,"UT",2) || !strncmp(VR,"UN",2)){
	
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	if(sdum) {
	  fprintf(stderr,"inconsistency in VR field, should read 0, instead is %0.4X\n",sdum);
	  exit(1);
	}
      }
      if(!strncmp(VR,"OB",2) || !strncmp(VR,"OW",2) || !strncmp(VR,"OF",2)
	 || !strncmp(VR,"SQ",2) || !strncmp(VR,"UT",2) || !strncmp(VR,"UN",2)
	 || !VR[0]){
	
	if(fread(&length,4,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
      }
      else {
	if(fread(&slength,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	length = (long)slength;
      }
      
      /* new sequence */

      if(!strncmp(VR,"SQ",2)) {
	if(verbose) {
	  indent();
	  printf("tag: group %0.4X, element %0.4X VR %c%c length %ld (%0.8lX)\n",mytag.group,mytag.element,VR[0],VR[1],length,length);
	}
	sqstart = ftell(*fpp);
	sqlength = length;
	for(item = 0;;item++) {
	  if((sqlength != -1) && (ftell(*fpp) - sqstart) >= sqlength) {
	    if(verbose)
	      printf("end of SQ list after %d items.  %ld bytes read, sqlength = %ld\n",item,ftell(*fpp)-sqstart,sqlength);
	    break;
	  }
	  if(fread(&slength,2,1,*fpp) != 1) {
	    perror("read error");
	    exit(1);
	  }
	  if(slength != 0xFFFE) {
	    indent();
	    fprintf(stderr,"malformed list after SQ VR, expected 0xFFFE, got 0x%X\n",slength);
	    exit(1);
	  }
	  if(fread(&slength,2,1,*fpp) != 1) {
	    perror("read error");
	    exit(1);
	  }
	  if(slength == 0xE0DD) {
	    if(verbose) {
	      indent();
	      printf("end of SQ list after %d items.\n",item);
	    }
	    if(fread(&length,4,1,*fpp) != 1) {
	      perror("read error");
	      exit(1);
	    }
	    if(length) {
	      if(verbose) {
		indent();
		printf("unexpected length %ld after sequence delimiter\n",length);
	      }
	    }
	    if(fseek(*fpp,length,SEEK_CUR) == -1) {
	      indent();
	      printf("fseek failed...\n");
	      exit(1);
	    }
	    break;
	  }
	  else if(slength == 0xE000) {
	    if(verbose) {
	      indent();
	      printf("start of item %d\n",item);
	    }
	    if(fread(&length,4,1,*fpp) != 1) {
	      perror("read error");
	      exit(1);
	    }
	    level++;
	    iswantedloop(fpp,length);
	    level--;
	  }
	  else {
	    indent();
	    printf("unexpected item tag 0xFFFE, 0x%.X\n",slength);
	    exit(1);
	  }
	}
	continue;
      }

      switch(mytag.group) {
      case 0x20:
	switch(mytag.element) {
	case 0x11:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  thisseries = atoi(cbuf);
	  for(i=0; i<nseries; i++)
	    if(thisseries == series[i]) {
	      iwantit++;
	      break;
	    }
	  iswantedloopdone++;
	  break;
	default:
	  if(fseek(*fpp,length,SEEK_CUR) == -1) {
	    perror("fseek error 1");
	    exit(1);
	  }
	  break;
	}
	break;
      case 0x7fe0:
	switch(mytag.element) {
	case 0x10:
	  fprintf(stderr,"The sky is falling!  Reached end of dicom file without getting a series number!\n");
	  iswantedloopdone++;
	  continue;
	default:
	  if(fseek(*fpp,length,SEEK_CUR) == -1) {
	    perror("fseek error 2");
	    exit(1);
	  }
	  break;
	}
	break;
      case 0xFFFE:
	switch(mytag.element) {
	case 0xE000:  // beginning of a new data sequence
	  level++;
	  iswantedloop(fpp,length);
	  level--;
	  break;
	case 0xE00D:  // end of item
	  if(length) {
	    indent();
	    fprintf(stderr,"non-zero length for end-of-item, bailing...\n");
	    exit(1);
	  }
	  if(level)
	    return;
	  else
	    if(verbose) {
	      indent();
	      printf("got end-of-item at level 0...\n");
	    }
	  break;
	case 0xE0DD:  // end of sequence, shouldn't get this
	  indent();
	  fprintf(stderr,"got end-of-sequence without beginning of sequence, bailing...\n");
	  exit(1);
	  break; 
	default:
	  if(fseek(*fpp,length,SEEK_CUR) == -1) {
	    perror("fseek error");
	    exit(1);
	  }
	  break;
	}
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  printf("length = %d\n",length);
	  perror("fseek error 3");
	  exit(1);
	}
	break;
	
      }   
      if(itemlength != -1) {  // -1 means no length specified, going until end-of-item
	if((ftell(*fpp) - start) >= itemlength) {
	  if(verbose) {
	    indent();
	    printf("returning after reading %ld bytes, itemlength = %ld\n",ftell(*fpp)-start,itemlength);
	  }
	  return;
	}
      }
  }
}

readdicom(ip,fp)
     struct image *ip;
     FILE *fp;
{
  struct tag mytag;
  short sdum,slength;
  long length;
  char cbuf[256];

  fseek(fp,132L,SEEK_SET);
  dicomreadloopdone = 0;
  readdicomloop(&ip,&fp,-1);


}
readdicomloop(ipp,fpp,itemlength)
     struct image **ipp;
     FILE **fpp;
     long itemlength;
{
  long start;
  long sqstart,sqlength;
  struct tag {
    unsigned short group;
    unsigned short element;
  } mytag;
  char VR[2];
  unsigned int length;
  unsigned short slength;
  unsigned short sdum;
  unsigned int intdum;
  char cbuf[256];
  char *cptr;
  int item;
  char typestring[256];

  if(itemlength == 0)
    return;
  start = ftell(*fpp);
  for(; !dicomreadloopdone;) {
    if(fread(&mytag,sizeof(struct tag),1,*fpp) != 1) {
      if(feof(*fpp))
	return;
      perror("read error");
      exit(1);
    }
    if(verbose) {
      indent();
      printf("group %0.4X, element %0.4X\n",mytag.group,mytag.element);
    }
    if(mytag.group == 0xFFFE)
      VR[0] = VR[1] = '\0';
    else
      if(fread(&VR[0],sizeof(VR),1,*fpp) != 1) {
	if(feof(*fpp))
	  return;
	perror("read error");
	exit(1);
      }
    
    if(!strncmp(VR,"OB",2) || !strncmp(VR,"OW",2) || !strncmp(VR,"OF",2)
       || !strncmp(VR,"SQ",2) || !strncmp(VR,"UT",2) || !strncmp(VR,"UN",2)){
      
      if(fread(&sdum,2,1,*fpp) != 1) {
	perror("read error");
	exit(1);
      }
      if(sdum) {
	fprintf(stderr,"inconsistency in VR field, should read 0, instead is %0.4X\n",sdum);
	exit(1);
      }
    }
    if(!strncmp(VR,"OB",2) || !strncmp(VR,"OW",2) || !strncmp(VR,"OF",2)
       || !strncmp(VR,"SQ",2) || !strncmp(VR,"UT",2) || !strncmp(VR,"UN",2)
       || !VR[0]){
      
      if(fread(&length,4,1,*fpp) != 1) {
	perror("read error");
	exit(1);
      }
    }
    else {
      if(fread(&slength,2,1,*fpp) != 1) {
	perror("read error");
	exit(1);
      }
      length = (long)slength;
    }
    
    /* new sequence */
    
    if(!strncmp(VR,"SQ",2)) {
      if(verbose) {
	indent();
	printf("tag: group %0.4X, element %0.4X VR %c%c length %ld (%0.8lX)\n",mytag.group,mytag.element,VR[0],VR[1],length,length);
      }
      sqstart = ftell(*fpp);
      sqlength = length;
      for(item = 0;;item++) {
	if((sqlength != -1) && (ftell(*fpp) - sqstart) >= sqlength) {
	  if(verbose)
	    printf("end of SQ list after %d items.  %ld bytes read, sqlength = %ld\n",item,ftell(*fpp)-sqstart,sqlength);
	  break;
	}
	if(fread(&slength,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	if(slength != 0xFFFE) {
	  indent();
	  fprintf(stderr,"malformed list after SQ VR, expected 0xFFFE, got 0x%X\n",slength);
	  exit(1);
	}
	if(fread(&slength,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	if(slength == 0xE0DD) {
	  if(verbose) {
	    indent();
	    printf("end of SQ list after %d items.\n",item);
	  }
	  if(fread(&length,4,1,*fpp) != 1) {
	    perror("read error");
	    exit(1);
	  }
	  if(length) {
	    if(verbose) {
	      indent();
	      printf("unexpected length %ld after sequence delimiter\n",length);
	    }
	  }
	  if(fseek(*fpp,length,SEEK_CUR) == -1) {
	    indent();
	    printf("fseek failed...\n");
	    exit(1);
	  }
	  break;
	}
	else if(slength == 0xE000) {
	  if(verbose) {
	    indent();
	    printf("start of item %d\n",item);
	  }
	  if(fread(&length,4,1,*fpp) != 1) {
	    perror("read error");
	    exit(1);
	  }
	  level++;
	  readdicomloop(ipp,fpp,length);
	  level--;
	}
	else {
	  indent();
	  printf("unexpected item tag 0xFFFE, 0x%.X\n",slength);
	  exit(1);
	}
      }
      continue;
    }
    switch(mytag.group) {
    case 0x2:
      switch(mytag.element) {
      case 0x3:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->mediasopUID = malloc(length+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->mediasopUID,cbuf);
	break;
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error 4");
	  exit(1);
	}
	break;
      }
      break;
    case 0x8:
      switch(mytag.element) {
      case 0x1:
	printf("tag \"length to end\" has length %ld\n",length);
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	printf("cbuf = %s\n",cbuf);
	break;
      case 0x8:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	strcpy(typestring,cbuf);
	strtok(cbuf,"\\");
	strtok((char *)0,"\\");
	cptr = strtok((char *)0,"\\");
	if((cptr == (char *)0)) {
	  fprintf(stderr,"image type string has unknown format: %s\n",typestring);
	  (*ipp)->imtype = 'X';
	}
	else
	  (*ipp)->imtype = cptr[0];
	break;
      case 0x18:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->sopUID = malloc(length+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->sopUID,cbuf);
	break;
      case 0x20:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->date = malloc(length+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->date,cbuf);
	if((*ipp)->date[strlen((*ipp)->date)-1] == ' ')  // remove possible trailing space
	  (*ipp)->date[strlen((*ipp)->date)-1] = '\0';
	break;
      case 0x32:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->acqtime = malloc(length + 1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->acqtime,cbuf);
	break;
      case 0x50:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->accession = atol(cbuf);
	break;
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error 5");
	  exit(1);
	}
	break;
      }
      break;
    case 0x10:
      switch(mytag.element) {
      case 0x10:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->name = malloc(length+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->name,cbuf);
	if((*ipp)->name[strlen((*ipp)->name)-1] == ' ')
	  (*ipp)->name[strlen((*ipp)->name)-1] = '\0';
	break;
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error");
	  exit(1);
	}
	break;
      }
      break;
    case 0x18:
      switch(mytag.element) {
      case 0x50:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->slicethickness = malloc(length + 1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->slicethickness,cbuf);
	break;
      case 0x80:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->tr = atof(cbuf);
	break;
      case 0x81:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->te = atof(cbuf);
	break;
      case 0x84:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->larmor_freq = atof(cbuf);
	break;
      case 0x86:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->echono = atoi(cbuf);
	break;
      case 0x1314:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->flipangle = (float)atof(cbuf);
	break;
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error");
	  exit(1);
	}
	break;
      }
      break;
    case 0x20:
      switch(mytag.element) {
      case 0xe:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->seriesUID = malloc(length+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->seriesUID,cbuf);
	break;
      case 0x11:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->series = atoi(cbuf);
	break;
      case 0x13:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	(*ipp)->imageno = atoi(cbuf);
	break;
      case 0x32:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->imagepos = malloc(length + 1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->imagepos,cbuf);
	break;
      case 0x1041:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->sliceloc = malloc(length + 1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->sliceloc,cbuf);
	if((*ipp)->sliceloc[strlen((*ipp)->sliceloc)-1] == ' ')
	  // yet another goofiness of the DICOM format
	  (*ipp)->sliceloc[strlen((*ipp)->sliceloc)-1] = '\0';
	break;
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error");
	  exit(1);
	}
	break;
      }
      break;
    case 0x28:
      switch(mytag.element) {
      case 0x10:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	if(!level)  // nested objects have rows attribute with same tag!
	  (*ipp)->rows = (int)sdum;
	break;
      case 0x11:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	if(!level)  // nested objects have cols attribute with same tag!
	  (*ipp)->cols = (int)sdum;
	break;
      case 0x30:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	if(((*ipp)->pixelspacing = malloc(length + 1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy((*ipp)->pixelspacing,cbuf);
	break;
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error");
	  exit(1);
	}
	break;
      }
      break;
    case 0x7fe0:
      switch(mytag.element) {
      case 0x10:
	if(length == -1) {  // list of encapsulated image data
	  for(;;) {
	    if(fread(&mytag,sizeof(struct tag),1,*fpp) != 1) {
	      if(feof(*fpp))
		return;
	      perror("read error");
	      exit(1);
	    }
	    if((mytag.group == 0xFFFE) && (mytag.element == 0xE0DD)) { // end of list
	      if(verbose) {
		indent();
		printf("end of encapsulated image data\n");
	      }
	      break;
	    }
	    if((mytag.group != 0xFFFE) || (mytag.element != 0xE000)) {
	      fprintf(stderr,"malformed encapsulated image data list, bailing...\n");
	      fprintf(stderr,"next tag: group %x, element %x\n",mytag.group,mytag.element);
	      exit(1);
	    }
	    if(fread(&length,4,1,*fpp) != 1) {
	      perror("read error");
	      exit(1);
	    }
	    if(verbose) {
	      indent();
	      printf("encapsulated image data list, length = %d\n",length);
	    }
	    if(fseek(*fpp,length,SEEK_CUR) == -1) {
	      perror("fseek error");
	      exit(1);
	    }
	  }
	  break;
	}
	else {
	  if(!level) {
	    if(((*ipp)->image = (short *)malloc(length)) == NULL) {
	      perror("malloc error");
	      exit(1);
	    }
	    if(fread((*ipp)->image,1,length,*fpp) != length) {
	      perror("read error in picture data");
	      exit(1);
	    }
	    dicomreadloopdone = 1;
	  } else
	    if(fseek(*fpp,length,SEEK_CUR) == -1) {
	      perror("fseek error");
	      exit(1);
	    }
	  break;
	}
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error");
	  exit(1);
	}
	break;
      }
      break;
    case 0xFFFE:
      switch(mytag.element) {
      case 0xE000:  // beginning of a new data sequence
	level++;
	readdicomloop(ipp,fpp,length);
	level--;
	break;
      case 0xE00D:  // end of item
	if(length) {
	  indent();
	  fprintf(stderr,"non-zero length for end-of-item, bailing...\n");
	  exit(1);
	}
	if(level)
	  return;
	else
	  if(verbose) {
	    indent();
	    printf("got end-of-item at level 0...\n");
	  }
	break;
      case 0xE0DD:  // end of sequence, shouldn't get this
	indent();
	fprintf(stderr,"got end-of-sequence without beginning of sequence, bailing...\n");
	exit(1);
	break; 
      default:
	if(fseek(*fpp,length,SEEK_CUR) == -1) {
	  perror("fseek error");
	  exit(1);
	}
	break;
      }
    default:
      if(fseek(*fpp,length,SEEK_CUR) == -1) {
	perror("fseek error");
	exit(1);
      }
      break;
    }
    if(itemlength != -1) {  // -1 means no length specified, going until end-of-item
      if((ftell(*fpp) - start) >= itemlength) {
	if(verbose) {
	  indent();
	  printf("returning after reading %ld bytes, itemlength = %ld\n",ftell(*fpp)-start,itemlength);
	}
	return;
      }
    }
  }
}


getseries(argv)
     char *argv;
{
  int i,j,k,ncommas;
  char *endp;

  ncommas = 0;
  for(i=0;argv[i] != '\0';i++) 
    if(argv[i] == ',')
      ncommas++;
  nseries = ncommas + 1;
  if((series = (int *)malloc(sizeof(int) * nseries)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0,j=0,k=0; k < nseries ;i++) {
    if((argv[i] == '\0') || (argv[i] == ',')) {
      if(j==i) {
	usage((char **)0);
	exit(1);
      }
      argv[i] = '\0';
      series[k] = strtol(&argv[j],&endp,10);
      if(endp != (&argv[i])) {
	fprintf(stderr,"screwy series string went bad at \"%s\"\n",endp);
	usage((char **)0);
      }
      j = i + 1;
      k++;
    }
  }
}
