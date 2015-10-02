#include<stdlib.h>
#include<stdio.h>
#include<errno.h>
#include<dirent.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int32_t bydate, copyonly, npatnames, ndates;
char **patname,**date;
char destdir[PATH_MAX+1];
struct dirnames {
  char *name;
  int32_t nfiles;
  int32_t created;
  struct dirnames *next;
} basedirname;

int32_t level,verbose;
char *cp,*datep;
int32_t series,gotprot,gotpatname,gotdate;

main(argc,argv)
     int32_t argc;
     char **argv;
{
  int32_t i;
  struct dirnames *dnp;
  char *meebs;

  if(argc == 1) {
    fprintf(stderr,"Program for separating DICOM files into directories by patient name\n");
    fprintf(stderr,"usage: %s [-d destdir] [-bydate] dir1 [dir2 .. ]\n",argv[0]);
    fprintf(stderr,"\t -d destdir \troot directory to put patient name subdirectories in\n");
    fprintf(stderr,"\tdir1, dir2, etc are directories containing dicom files.\n");
    fprintf(stderr,"\t\tThese will be searched recursively.\n");
    fprintf(stderr,"\t-bydate\tmake separate subdirectories by study date.\n");
    fprintf(stderr,"\t-c\tcopy files to directory but don't try to delete originals.\n");
    fprintf(stderr,"\t-v\tverbose\n");
    exit(1);
  }
  /*  printf("DT_UNKNOWN: %xd\n",DT_UNKNOWN);
  printf("DT_REG: %xd\n",DT_DIR);
  printf("DT_DIR: %xd\n",DT_DIR);
  printf("DT_FIFO: %xd\n",DT_FIFO);
  printf("DT_SOCK: %xd\n",DT_SOCK);
  printf("DT_CHR: %xd\n",DT_CHR);
  printf("DT_BLK: %xd\n",DT_BLK);
  */

  level = 0;
  npatnames = ndates = bydate = copyonly = 0;
  patname = date = (char **)0;
  basedirname.name = (char *)0;
  destdir[0] = '\0';
  for(i=1; i < argc; i++) {
    if(!strcmp(argv[i],"-bydate")) {
      bydate++;
      continue;
    }
    if(!strcmp(argv[i],"-d")) {
      i++;
      strcpy(destdir,argv[i]);
      printf("destination directory: %s\n",destdir);
      continue;
    }
    if(!strcmp(argv[i],"-c")) {
      copyonly++;
      continue;
    }
    if(!strcmp(argv[i],"-v")) {
      verbose++;
      continue;
    }
    dodirectory(argv[i]);
  }
  if(copyonly)
    meebs = "Copied";
  else
    meebs = "Moved";

    for(dnp = &basedirname; dnp; dnp = dnp->next) {
      if(dnp->created)
	printf("%s %d files to new directory %s\n",meebs,dnp->nfiles,dnp->name);
      else
	printf("%s %d files to existing directory %s\n",meebs,dnp->nfiles,dnp->name);
    }
}


dodirectory(dirname)
char *dirname;
{
  struct stat sbuf;
  FILE *fp;
  char fourchar[5];
  struct tag {
    int16_t group;
    int16_t element;
    char VR[2];
  } mytag;
  int16_t sdum,slength;
  int32_t length;
  char cbuf[256];
  char sourcefile[PATH_MAX+1];
  char destfile[PATH_MAX+1];
  char dtemp[PATH_MAX+1];
  char prot[256];
  char sysstring[3*PATH_MAX];
  int32_t first;
  struct dirent **direntp,*dp;
  int32_t j;
  int32_t nents;
  struct dirnames *dnp;

  printf("directory %s:\n",dirname);
  fflush(stdout);
  if((nents = scandir(dirname,&direntp,0,0)) < 0) {
    perror("scandir failed");
    fprintf(stderr,"directory name = %s\n",dirname);
    exit(1);
  }
  for(j=0,dp = direntp[0];j < nents;dp = direntp[++j]) {
    cp = datep = (char *)0;
    gotpatname = gotdate = 0;
    
    // stat() knows that iso9660 directories are directories.  scandir calls both files
    // and directories on cd rom type "unknown."  Must use stat() to see if this is a directory.
    
    sprintf(sourcefile,"%s/%s",dirname,dp->d_name);
    if(stat(sourcefile,&sbuf)) {
      fprintf(stderr,"stat failed for file %s\n",sourcefile);
      perror("");
      exit(1);
    }
    if(S_ISDIR(sbuf.st_mode)) {
      if(!strcmp(dp->d_name,".") || !strcmp(dp->d_name,".."))
	continue;
      dodirectory(sourcefile);
      continue;
    }
    // files on CDROM seem to be type "unknown"
    //    if((dp->d_type != DT_REG) && (dp->d_type != DT_UNKNOWN))
    /*    printf("File type for file %s:                ",sourcefile);
           switch (sbuf.st_mode & S_IFMT) {
           case S_IFBLK:  printf("block device\n");            break;
           case S_IFCHR:  printf("character device\n");        break;
           case S_IFDIR:  printf("directory\n");               break;
           case S_IFIFO:  printf("FIFO/pipe\n");               break;
           case S_IFLNK:  printf("symlink\n");                 break;
           case S_IFREG:  printf("regular file\n");            break;
           case S_IFSOCK: printf("socket\n");                  break;
           default:       printf("unknown?\n");                break;
           }
    */
   if(!S_ISREG(sbuf.st_mode)) {
      printf("skipping non-regular file %s\n",sourcefile);
      continue;
    }
    if((fp = fopen(sourcefile,"r")) == NULL) {
      fprintf(stderr,"could not open file %s: ",sourcefile);
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
    if(strcmp(fourchar,"DICM")) {  // not DICOM file
      continue;
    }
    gotdate = gotpatname = 0;
    loop(&fp,-1);
    fclose(fp);
    if(!cp) {
      printf("no patient name in file %s, did not sort\n",dp->d_name);
      continue;
    }
    if(bydate && (!datep)) {
      printf("no date field in file %s, did not sort\n",dp->d_name);
      continue;
    }
    
    // add to list of directories and increment number of files
    
    if(destdir[0])
      if(bydate)
	sprintf(destfile,"%s/%s/%s\0",destdir,cp,datep);
      else
	sprintf(destfile,"%s/%s\0",destdir,cp);
    else
      if(bydate)
	sprintf(destfile,"./%s/%s\0",cp,datep);
      else
	sprintf(destfile,"./%s/0",cp);
    if(basedirname.name == (char *)0) {
      if((basedirname.name = malloc(strlen(destfile)+1)) == NULL) {
	perror("malloc failed");
	exit(1);
      }
      strcpy(basedirname.name,destfile);
      basedirname.nfiles = 0;
      basedirname.next = (struct dirnames *)0;
      basedirname.created = 0;
    }
    for(dnp = &basedirname; dnp && strcmp(destfile,dnp->name); dnp = dnp->next);
    if(!dnp) {
      for(dnp = &basedirname; dnp->next; dnp = dnp->next);
      if((dnp->next = (struct dirnames *)malloc(sizeof(struct dirnames))) == NULL) {
	perror("malloc failed");
	exit(1);
      }
      dnp = dnp->next;
      dnp->next = (struct dirnames *)0;
      if((dnp->name = malloc(strlen(destfile)+1)) == NULL) {
	perror("malloc failed");
	exit(1);
      }
      strcpy(dnp->name,destfile);
      dnp->nfiles = 0;
      dnp->created = 0;
    }
    
    // create directory if necessary

    if(destdir[0])
      sprintf(destfile,"%s/%s",destdir,cp);
    else
      sprintf(destfile,"./%s",cp);
    if(stat(destfile,&sbuf)) {
      printf("making directory %s\n",destfile);
      if(mkdir(destfile,0777)) {
	fprintf(stderr,"mkdir failed for directory %s: ",destfile);
	perror("");
	exit(1);
      }
      dnp->created++;
    }
    else if(!S_ISDIR(sbuf.st_mode)) {
      fprintf(stderr,"error: %s exists and is not a directory: ",destfile);
      perror("");
      exit(1);
    }
    if(bydate) {
      if(destdir[0])
	sprintf(destfile,"%s/%s/%s",destdir,cp,datep);
      else
	sprintf(destfile,"./%s/%s",cp,datep);
      if(stat(destfile,&sbuf)) {
	printf("making directory %s\n",destfile);
	if(mkdir(destfile,0777)) {
	  fprintf(stderr,"mkdir failed for directory %s: ",destfile);
	  perror("");
	  exit(1);
	}
	dnp->created++;
      }
      else if(!S_ISDIR(sbuf.st_mode)) {
	fprintf(stderr,"error: %s exists and is not a directory: ",destfile);
	perror("");
	exit(1);
      }
    }
    
    if(destdir[0])
      if(bydate)
	sprintf(destfile,"%s/%s/%s/%s\0",destdir,cp,datep,dp->d_name);
      else
	sprintf(destfile,"%s/%s/%s\0",destdir,cp,dp->d_name);
    else
      if(bydate)
	sprintf(destfile,"./%s/%s/%s\0",cp,datep,dp->d_name);
      else
	sprintf(destfile,"./%s/%s\0",cp,dp->d_name);
    if(!stat(destfile,&sbuf))
      printf("file %s exists, skipping...\n",destfile);
    else {
      if(copyonly)
	sprintf(sysstring,"cp -f \"%s\" \"%s\"\0",sourcefile,destfile);
	else
	  sprintf(sysstring,"mv -f \"%s\" \"%s\"\0",sourcefile,destfile);
      if(system(sysstring)) {   // using system because rename() doesn't work across filesystems
	perror("system call failed");
	exit(1);
      }
      dnp->nfiles++;
    }
  }
}

loop(fpp,itemlength)
     FILE **fpp;
     int32_t itemlength;
{
  int32_t start;
  int32_t sqstart,sqlength;
  int32_t done,curgroup;
  struct tag {
    u_int16_t group;
    u_int16_t element;
  } mytag;
  char VR[2];
  u_int32_t length;
  u_int16_t slength;
  u_int16_t sdum;
  u_int32_t intdum;
  char cbuf[256];
  int32_t item;
  int32_t picnum;
  char picfile[256];
  int32_t k;

  if(itemlength == 0)
    return;
  start = ftell(*fpp);
  curgroup = -1;
  for(;(!gotdate) || (!gotpatname);) {
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
	length = (int32_t)slength;
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
	    if(verbose) {
	      indent();
	      printf("end of SQ list after %d items.  %ld bytes read, sqlength = %ld\n",item,ftell(*fpp)-sqstart,sqlength);
	    }
	    if(itemlength != -1) {
	      if((ftell(*fpp) - start) >= itemlength)
		return;
	    }
	    break;
	  }
	  if(fread(&slength,2,1,*fpp) != 1) {
	    perror("read error");
	    exit(1);
	  }
	  if(slength != 0xFFFE) {
	    indent();
	    printf("malformed list after SQ VR, expected 0xFFFE, got 0x%X\n",slength);
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
	      indent();
	      printf("unexpected length %ld after sequence delimiter\n",length);
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
	    loop(fpp,length);
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

      if(verbose) {
	indent();
	printf("tag: group %0.4X, element %0.4X VR %c%c length %ld (%0.8lX)\n",mytag.group,mytag.element,VR[0],VR[1],length,length);
      }
      switch(mytag.group) {
      case 0x8:
	switch(mytag.element) {
	case 0x20:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  if(cbuf[length-1] == ' ')  // strip away trailing space
	    cbuf[length-1] = '\0';
	  for(k=0; k<ndates; k++) {
	    if(!strcmp(date[k],cbuf)) {
	      datep = date[k];
	      break;
	    }
	  }
	  if(k == ndates) {
	    ndates++;
	    if(date == (char **)0) {
	      if((date = (char **)malloc(sizeof(char *))) == NULL) {
		perror("malloc failed");
		exit(1);
	      }
	    }
	    else
	      if((date = (char **)realloc(date,sizeof(char *)*ndates)) == NULL) {
		perror("realloc failed");
		exit(1);
	      }
	    if((date[ndates - 1] = malloc(strlen(cbuf)+1)) == NULL) {
	      perror("malloc failed");
	      exit(1);
	    }
	    strcpy(date[ndates - 1],cbuf);
	    datep = date[ndates - 1];
	  }
	  gotdate++;
	  break;
	default:
	  if(fseek(*fpp,length,SEEK_CUR) == -1) {
	    perror("fseek error");
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
	  if(cbuf[length-1] == ' ')  // strip away trailing space
	    cbuf[length-1] = '\0';
	  for(k=0; k<npatnames; k++) {
	    if(!strcmp(patname[k],cbuf)) {
	      cp = patname[k];
	      break;
	    }
	  }
	  if(k == npatnames) {
	    npatnames++;
	    if(patname == (char **)0) {
	      if((patname = (char **)malloc(sizeof(char *))) == NULL) {
		perror("malloc failed");
		exit(1);
	      }
	    }
	    else
	      if((patname = (char **)realloc(patname,sizeof(char *)*npatnames)) == NULL) {
		perror("realloc failed");
		exit(1);
	      }
	    if((patname[npatnames - 1] = malloc(strlen(cbuf)+1)) == NULL) {
	      perror("malloc failed");
	      exit(1);
	    }
	    strcpy(patname[npatnames - 1],cbuf);
	    cp = patname[npatnames - 1];
	  }
	  gotpatname++;
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
	    for(picnum = 0;;picnum++) {
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
		printf("next tag: group %x, element %x\n",mytag.group,mytag.element);
		exit(1);
	      }
	      if(fread(&length,4,1,*fpp) != 1) {
		perror("read error");
		exit(1);
	      }
	      if(fseek(*fpp,length,SEEK_CUR) == -1) {
		perror("fseek error");
		exit(1);
	      }
	      if(verbose) {
		indent();
		printf("encapsulated image data list, length = %d\n",length);
	      }
	    }
	    break;
	  }
	  else {
	    if(fseek(*fpp,length,SEEK_CUR) == -1) {
	      perror("fseek error");
	      exit(1);
	    }
	    if(verbose) {
	      indent();
	      printf("non-encapsulated image, length = %d\n",length);
	    }
	  }
	  break;
	  
	default:
	  break;
	}
	break;
      case 0xFFFE:
	switch(mytag.element) {
	case 0xE000:  // beginning of a new data sequence
	  level++;
	  loop(fpp,length);
	  level--;
	  break;
	case 0xE00D:  // end of item
	  if(length) {
	    indent();
	    printf("non-zero length for end-of-item, bailing...\n");
	    exit(1);
	  }
	  if(level)
	    return;
	  else
	    indent();
	    printf("got end-of-item at level 0...\n");
	  break;
	case 0xE0DD:  // end of sequence, shouldn't get this
	  indent();
	  printf("got end-of-sequence without beginning of sequence, bailing...\n");
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

indent()
{
  int32_t i;

  for(i=0; i<level; i++)
    printf("\t");
}
