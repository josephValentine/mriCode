/* more sophisticated DICOM dump program, deals with SQ lists.
 */

#include<stdlib.h>
#include<stdio.h>
#include<errno.h>

int32_t level,verbose;
int16_t *pixbuf;
FILE *ofp;

main(argc,argv)
     int32_t argc;
     char **argv;

{

  int32_t file;
  FILE *fp;

  if(argc == 1) {
    fprintf(stderr,"usage: %s [-v] file1 [file2 file3 ...]\n",argv[0]);
    fprintf(stderr,"\t-v\tverbose mode, reports all dicom tags\n");
    fprintf(stderr,"\t-v\tadditional -v's mean more verbose\n");
    fprintf(stderr,"\tfilen\tlist of files to dump\n");
    exit(1);
  }

  verbose = 0;
  for(file=1; file<argc; file++) {
    if(!strcmp(argv[file],"-v")) {
      verbose++;
      continue;
    }
    init(&fp,argv[file]);
    loop(&fp,-1);
  }
}

init(fpp,fname)
     FILE **fpp;
     char *fname;
{
  char fourchar[5];

  if(((*fpp) = fopen(fname,"r")) == NULL) {
    fprintf(stderr,"can't open %s for reading: ",fname);
      perror("");
      exit(1);
  }
  printf("file %s:\n",fname);
  // This section was commented out, can't remember why.  Some DICOM noncompliance issue.
  if(fseek(*fpp,128L,SEEK_CUR) == -1) {
    perror("fseek error");
    exit(1);
  }
  if(fread(fourchar,1,4,*fpp) != 4) {
    perror("read error");
    exit(1);
  }
  fourchar[4] = '\0';
  if(strcmp(fourchar,"DICM")) {
    fprintf(stderr,"%s does not appear to be a DICOM file\n",fname);
    // exit(0);
  }
  // end of commented out section.  Following line was in.
  /*    fseek(fp,16L,SEEK_SET); */
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

  if(itemlength == 0)
    return;
  start = ftell(*fpp);
  done = 0;
  curgroup = -1;
  for(;!done;) {
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
	  if(verbose > 1) {
	    indent();
	    printf("currently at length %ld in this SQ list...\n",ftell(*fpp)-sqstart);
	  }
	  if((sqlength != -1) && (ftell(*fpp) - sqstart) >= sqlength) {
	    if(verbose) {
	      indent();
	      printf("end of SQ list after %d items.  %ld bytes read, sqlength = %ld\n",item,ftell(*fpp)-sqstart,sqlength);
	    }
	    if(verbose > 1) {
	      indent();
	      printf("currently at length %ld in this item...\n",ftell(*fpp) - start);
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
	    indent();
	    printf("end of SQ list after %d items.\n",item);
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
	    if(fread(&length,4,1,*fpp) != 1) {
	      perror("read error");
	      exit(1);
	    }
	    indent();
	    printf("start of item %d, length = %ld\n",item,length);
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
	indent();
	printf("Group %0.4X:\n",mytag.group);
	curgroup = mytag.group;
      }
      if(verbose) {
	indent();
	printf("tag: group %0.4X, element %0.4X VR %c%c length %ld (%0.8lX)\n",mytag.group,mytag.element,VR[0],VR[1],length,length);
      }
      switch(mytag.group) {
      case 0x2:
	switch(mytag.element) {
	case 0x0:
	  if(fread(&intdum,sizeof(intdum),1,*fpp) != 1) {
	    perror("read error");
	    exit(1);
	  }
	  indent();
	  printf("**Group length: %d\n",intdum);
	  break;
	case 0x1:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**File meta information version: %s\n",cbuf);
	  break;
	case 0x2:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**Media storage SOP class UID: %s\n",cbuf);
	  break;
	case 0x3:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**Media storage SOP instance UID: %s\n",cbuf);
	  break;
	case 0x10:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**Transfer syntax UID: %s\n",cbuf);
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
	  indent();
	  printf("**image type: %s\n",cbuf);
	  break;
	case 0x16:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**SOP class UID: %s\n",cbuf);
	  break;
	case 0x18:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**SOP instance UID: %s\n",cbuf);
	  break;
	case 0x20:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**+study date: %s\n",cbuf);
	  break;
	case 0x30:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**+study time: %s\n",cbuf);
	  break;
	case 0x31:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("series time: %s\n",cbuf);
	  break;
	case 0x32:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("acquisition time: %s\n",cbuf);
	  break;
	case 0x33:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("content time: %s\n",cbuf);
	  break;
	case 0x50:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("accession number: %s\n",cbuf);
	  break;
	case 0x70:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("manufacturer: %s\n",cbuf);
	  break;
 	case 0x1030:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("study description: %s\n",cbuf);
	  break;
	case 0x103e:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("series description: %s\n",cbuf);
	  break;
	case 0x9208:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**+complex image component: %s\n",cbuf);
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
	indent();
	printf("**+patient name: %s\n",cbuf);
	break;
      case 0x20:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+patient ID: %s\n",cbuf);
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
	case 0x20:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("scanning seqeunce: %s\n",cbuf);
	  break;
	case 0x21:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("scanning seqeunce modifier: %s\n",cbuf);
	  break;
	case 0x23:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("MR acquisition type: %s\n",cbuf);
	  break;
	case 0x24:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("sequence name: %s\n",cbuf);
	  break;
	case 0x50:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**+slice thickness: %s\n",cbuf);
	  break;
	case 0x80:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("TR: %s\n",cbuf);
	  break;
	case 0x81:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**+TE: %s\n",cbuf);
	  break;
	case 0x84:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("**+imaging frequency: %s\n",cbuf);
	  break;
	case 0x86:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("echo number(s): %s\n",cbuf);
	  break;
	case 0x1030:
	  if(fread(cbuf,1,length,*fpp) != length) {
	    perror("read error");
	    exit(1);
	  }
	  cbuf[length] = '\0';
	  indent();
	  printf("protocol name: %s\n",cbuf);
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
      case 0xD:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**study instance UID: %s...\n",cbuf);
	break;
      case 0xE:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**series instance UID: %s...\n",cbuf);
	break;
      case 0x10:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+series description: %s...\n",cbuf);
	break;
      case 0x11:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+series number: %s...\n",cbuf);
	break;
      case 0x13:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+instance number: %s...\n",cbuf);
	break;
      case 0x12:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+acquisition number: %s\n",cbuf);
	break;
      case 0x32:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+image position: %s\n",cbuf);
	break;
      case 0x37:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**image orientation: %s\n",cbuf);
	break;
      case 0x52:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**frame of reference UID: %s\n",cbuf);
	break;
      case 0x1041:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**+slice location: %s\n",cbuf);
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
      case 0x2:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**samples per pixel: %d\n",sdum);
	break;
      case 0x4:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**photometric interpretation: %s\n",cbuf);
	break;
      case 0x10:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**rows: %d\n",sdum);
	break;
      case 0x11:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**columns: %d\n",sdum);
	break;
      case 0x30:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("**pixel spacing: %s\n",cbuf);
	break;
      case 0x100:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**bits allocated: %d\n",sdum);
	break;
      case 0x101:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**bits stored: %d\n",sdum);
	break;
      case 0x102:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**high bit: %d\n",sdum);
	break;
      case 0x103:
	if(fread(&sdum,2,1,*fpp) != 1) {
	  perror("read error");
	  exit(1);
	}
	indent();
	printf("**pixel representation: %d\n",sdum);
	break;
      case 0x1050:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("window center: %s\n",cbuf);
	break;
      case 0x1051:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("window width: %s\n",cbuf);
	break;
      case 0x1052:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("rescale intercept: %s\n",cbuf);
	break;
      case 0x1053:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("rescale slope: %s\n",cbuf);
	break;
      case 0x1054:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("rescale type: %s\n",cbuf);
	break;
      case 0x1055:
	if(fread(cbuf,1,length,*fpp) != length) {
	  perror("read error");
	  exit(1);
	}
	cbuf[length] = '\0';
	indent();
	printf("window center and width explanation: %s\n",cbuf);
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
		indent();
		printf("end of encapsulated image data\n");
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
	      if((pixbuf = (int16_t *)malloc(length)) == NULL) {
		perror("malloc error");
		exit(1);
	      }
	      if(fread(pixbuf,1,length,*fpp) != length) {
		perror("read error in picture data");
		exit(1);
	      }
	      /*	      sprintf(picfile,"./pic%d\0",picnum);
	      if((ofp = fopen(picfile,"w")) == NULL) {
		perror("error opening file");
		exit(1);
	      }
	      if(fwrite(pixbuf,1,length,ofp) != length) {
		perror("write error");
		exit(1);
	      }
	      */
	      free(pixbuf);
	      indent();
	      printf("encapsulated image data list, length = %d\n",length);
	      /*	      if(fseek(*fpp,length,SEEK_CUR) == -1) {
		perror("fseek error");
		exit(1);
	      }
	      */
	    }
	    break;
	  }
	  else {
	    if((pixbuf = (int16_t *)malloc(length)) == NULL) {
	      perror("malloc error");
	      exit(1);
	    }
	    if(fread(pixbuf,1,length,*fpp) != length) {
	      perror("read error in picture data");
	      exit(1);
	    }
	    /*	    if((ofp = fopen("./pic","w")) == NULL) {
	      perror("error opening file");
	      exit(1);
	    }
	    if(fwrite(pixbuf,1,length,ofp) != length) {
	      perror("write error");
	      exit(1);
	    }
	    */
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
	if(verbose > 1) {
	  indent();
	  printf("currently at length %ld in this item...\n",ftell(*fpp) - start);
	}
	if((ftell(*fpp) - start) >= itemlength) {
	  indent();
	  printf("returning after reading %ld bytes, itemlength = %ld\n",ftell(*fpp)-start,itemlength);
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
