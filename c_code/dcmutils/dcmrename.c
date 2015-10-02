#include<stdlib.h>
#include<stdio.h>
#include<errno.h>
#include<dirent.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include<limits.h>

int32_t series;
int32_t image;
int32_t echo;
int32_t verbose,level,derived;
char imtype;
char *cptr;

main(argc,argv)
     int32_t argc;
     char **argv;
{
  int32_t i,j,k;
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
  int32_t echoname;
  int32_t phasename;
  struct stat sbuf;
  char destfile[PATH_MAX+1],sourcefile[PATH_MAX+1];

  if(argc == 1) {
    fprintf(stderr,"Program for renaming all DICOM files in directory by series and image number.\n");
    fprintf(stderr,"usage: %s [-e -p -v] dir1 [dir2 .. ]\n",argv[0]);
    fprintf(stderr,"\t-e\tinclude echo number in file name\n");
    fprintf(stderr,"\t-p\tseparate phase and magnitude images\n");
    fprintf(stderr,"\t-v\tverbose, includes some dicom header parsing output\n");
    exit(1);
  }
  echoname = phasename = 0;
  verbose = 0;
  level = 0;
  for(i=1; i < argc; i++) {
    if(!strcmp(argv[i],"-v")) {
      verbose++;
      continue;
    }
  if(!strcmp(argv[i],"-e")) {
      echoname = 1;
      continue;
    }
    if(!strcmp(argv[i],"-p")) {
      phasename = 1;
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
      derived = 0;
      sprintf(sourcefile,"%s/%s",argv[i],dp->d_name);
      if(verbose)
	printf("\nFile name: %s\n",sourcefile);
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
      image = -1;
      series = -1;
      echo = -1;
      imtype = '\0';
      loop(&fp,-1);

      fclose(fp);
      if(series == -1)
	printf("\nfailed to get series number for file %s\n\n",sourcefile);
      if(image == -1)
	printf("\nfailed to get image number for file %s\n\n",sourcefile);
      if(echo == -1)
	printf("\nfailed to get echo number for file %s\n\n",sourcefile);
      if(imtype == '\0')
	printf("\nfailed to get image type for file %s\n\n",sourcefile);
      sprintf(destfile,"%s/s%03di%03d",argv[i],series,image);
      if(echoname)
	sprintf(destfile,"%se%d",destfile,echo);
      if(phasename)
	sprintf(destfile,"%s%c",destfile,imtype);
      if(derived)
	sprintf(destfile,"%sd",destfile);
      if(!stat(destfile,&sbuf)) {
	if(strcmp(sourcefile,destfile))
	  printf("intended to rename %s to %s, but %s already exists.\n",sourcefile,destfile,destfile);
      }
      else 
	if(rename(sourcefile,destfile)) {
	  fprintf(stderr,"rename of %s to %s failed: ",sourcefile,destfile);
	  perror("");
	}
    }
  }
}

loop(fpp,itemlength)
     FILE **fpp;
     int32_t itemlength;
{
  int32_t start;
  int32_t sqstart,sqlength;
  int32_t curgroup;
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
  curgroup = -1;
  for(; (image == -1) || (series == -1) || (echo == -1) || (imtype == '\0');) {
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
      case 0x18:
	switch(mytag.element) {
	  case 0x86:
	    if(fread(cbuf,1,length,*fpp) != length) {
	      perror("read error");
	      exit(1);
	    }
	    cbuf[length] = '\0';
	    echo = atoi(cbuf);
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
	  cbuf[length] = '\0';
	  series = atoi(cbuf);
	  break;
	case 0x13:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  image = atoi(cbuf);
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
      case 0x8:
	switch(mytag.element) {
	case 0x8:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  if(!strncmp(cbuf,"DERIVED",7))
	    derived++;
	  strtok(cbuf,"\\");
	  strtok((char *)0,"\\");
	  cptr = strtok((char *)0,"\\");
	  if((!cptr) || (strlen(cptr) > 1)) {
	    fprintf(stderr,"image type string has unknown format: %s\n",cbuf);
	  }
	  //	  else 
	  imtype = cptr[0];
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
