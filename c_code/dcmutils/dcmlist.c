#include<stdlib.h>
#include<stdio.h>
#include<errno.h>
#include<dirent.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include<limits.h>

int32_t series,gotprot,gotdesc;
int32_t verbose,series_desc;
int32_t level;
char prot[256];
char desc[256];

main(argc,argv)
     int32_t argc;
     char **argv;
{
  int32_t i,j;
  struct dirent **direntp,*dp;
  FILE *fp;
  char fourchar[5];
  struct tag {
    int16_t group;
    int16_t element;
    char VR[2];
  } mytag;
  int32_t nents;
  int16_t sdum,slength;
  int32_t length;
  char cbuf[256];
  char sourcefile[PATH_MAX+1];
  struct bean {
    int32_t no;
    int32_t ims;
    char *prot;
    char *desc;
    struct bean *next;
    struct bean *prev;
  } *beanstar,*bp;
  int32_t minbean, lastbean;
  
  if(argc == 1) {
    fprintf(stderr,"Program for listing DICOM series numbers with protocol names.\n");
    fprintf(stderr,"usage: %s [-v] dir1 [dir2 .. ]\n",argv[0]);
    fprintf(stderr,"\t-v\tverbose\n");
    fprintf(stderr,"\t-s\tlist series description as well as protocol name.\n");
    exit(1);
  }
  beanstar = (struct bean *)0;
  verbose = series_desc = 0;
  for(i=1; i < argc; i++) {
    if(!strcmp(argv[i],"-v")) {
      verbose++;
      continue;
    }
    if(!strcmp(argv[i],"-s")) {
      series_desc++;
      continue;
    }
    if((nents = scandir(argv[i],&direntp,0,0)) < 0) {
      perror("scandir failed");
      exit(1);
    }
    printf("directory %s:\n",argv[i]);
    for(j=0,dp = direntp[j];dp;dp = direntp[++j]) {
      // files on CDROM seem to be type "unknown"
      if((dp->d_type != DT_REG) && (dp->d_type != DT_UNKNOWN))
	continue;
      sprintf(sourcefile,"%s/%s",argv[i],dp->d_name);
      if((fp = fopen(sourcefile,"r")) == NULL) {
	fprintf(stderr,"could not open file %s: ",sourcefile);
	perror("");
	exit(1);
      }
      if(verbose)
	printf("file %s\n",sourcefile);
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
	if(verbose)
	  printf("\tnot a dicom file\n");
	continue;
      }
      prot[0] = '\0';
      desc[0] = '\0';
      gotprot = gotdesc = 0;
      series = -1;
      level = 0;
      loop(&fp,-1);
      
      if(!beanstar) {
	if((beanstar = (struct bean *)malloc(sizeof(struct bean))) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	beanstar->prev = (struct bean *)0;
	beanstar->next = (struct bean *)0;
	beanstar->no = series;
	beanstar->ims = 0;
        beanstar->prot = (char *)0;
        beanstar->desc = (char *)0;
	bp = beanstar;
      }
      else {
	for(bp = beanstar; bp->next; bp = bp->next)
	  if(bp->no == series)
	    break;
	if(bp->no != series) {
	  if((bp->next = (struct bean *)malloc(sizeof(struct bean))) == NULL) {
	    perror("malloc failed");
	    exit(1);
	  }
	  bp->next->prev = bp;
	  bp->next->next = (struct bean *)0;
	  bp = bp->next;
	  bp->no = series;
	  bp->ims = 0;
	  bp->prot = (char *)0;
	  bp->desc = (char *)0;
	}
      }
      bp->ims++;
      if(bp->prot == (char *)0) {
	if((bp->prot = malloc(strlen(prot)+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
      }
      strcpy(bp->prot,prot);
      if(bp->desc == (char *)0) {
	if((bp->desc = malloc(strlen(desc)+1)) == NULL) {
	  perror("malloc failed");
	  exit(1);
	}
	strcpy(bp->desc,desc);
      }
      fclose(fp);
    }
    lastbean = -1;
    for(;;) {
      minbean = INT_MAX;
      for(bp = beanstar; bp; bp = bp->next) {
	if(bp->no > lastbean)
	  if(bp->no < minbean)
	    minbean = bp->no;
      }
      if(minbean == lastbean)
	break;
      if(minbean == INT_MAX)
	break;
      lastbean = minbean;
      for(bp = beanstar; (bp->no != minbean) && (bp); bp = bp->next);
      if(series_desc)
	printf("Series %d, %d images, %s (%s)\n",bp->no,bp->ims,bp->prot,bp->desc);
      else
	printf("Series %d, %d images, %s\n",bp->no,bp->ims,bp->prot);
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

  if(itemlength == 0)
    return;
  start = ftell(*fpp);
  done = 0;
  curgroup = -1;
  for(; (series == -1) || !gotprot || !gotdesc;) {
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
	    if(verbose)
	      printf("end of SQ list after %d items.  %ld bytes read, sqlength = %ld\n",item,ftell(*fpp)-sqstart,sqlength);
	    if(itemlength != -1) {  // -1 means no length specified, going until end-of-item
	      if((ftell(*fpp) - start) >= itemlength) {
		if(verbose) {
		  indent();
		  printf("returning after reading %ld bytes, itemlength = %ld\n",ftell(*fpp)-start,itemlength);
		}
		return;
	      }
	    }
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

      if(mytag.group != curgroup) {
	if(verbose) {
	  indent();
	  printf("\nGroup %0.4X:\n",mytag.group);
	}
	curgroup = mytag.group;
      }
      if(verbose) {
	indent();
	printf("tag: group %0.4X, element %0.4X VR %c%c length %ld (%0.8lX)\n",mytag.group,mytag.element,VR[0],VR[1],length,length);
      }
      switch(mytag.group) {
      case 0x08:
	switch(mytag.element) {
	case 0x103E:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  if(verbose) {
	    indent();
	    printf("Series description: %s\n",desc);
	  }
	  if(level == 0) {
	    cbuf[length] = '\0';
	    strcpy(desc,cbuf);
	    gotdesc++;
	  }
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
	case 0x1030:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  if(level == 0) {
	    cbuf[length] = '\0';
	    strcpy(prot,cbuf);
	    gotprot++;
	  }
	  if(verbose) {
	    indent();
	    printf("Protocol string: %s\n",prot);
	  }
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
      case 0x11:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	if(level == 0) {
	  cbuf[length] = '\0';
	  series = atoi(cbuf);
	}
	if(verbose) {
	  indent();
	  printf("Series number: %d\n",series);
	}
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
	    if(fseek(*fpp,length,SEEK_CUR) == -1) {
	      perror("fseek error");
	      exit(1);
	    }
	  }
	  break;
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
	  loop(fpp,length);
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


indent()
{
  int32_t i;

  for(i=0; i<level; i++)
    printf("\t");
}
