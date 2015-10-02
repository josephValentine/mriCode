// changed to create filename by image number rather than location
// changed to handle SQ lists, possibly nested, by just copying verbatim
// i.e. can't change anything in the SQ list


#include<stdlib.h>
#include<stdio.h>
#include<errno.h>
#include<string.h>
#include<limits.h>
#include"extern.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

// extern FILE *debug;

writeDicom(fname,ip,dbuf,seriesno,imageno,protstring,poff,cmdline)
     char *fname;
     struct image *ip;
     double *dbuf;
     int seriesno,imageno;
     char *protstring;
     double poff;
     char *cmdline;
{
  struct stat statbuf;
  struct tag tag;
  struct tag pofftag;
  short pofflen;
  short sdum,slength,*sbuf;
  short smax,smin;
  int i,done,oobs;
  int didmyfields;
  char poffstr[1024];
  int32_t length;
  int32_t cplength;
  char filename[1024];
  FILE *ofp,*ifp;
  union schar {
    short s;
    char c[sizeof(short)];
  } schar;
  char uid[256];
  char suid[256];
  float ran1();
  time_t dtime;

  char *cp = "It is rather for us to be here dedicated to the great task remaining before us -- that from these honored dead we take increased";
  
  if((sbuf = (short *)malloc(sizeof(short)*ip->rows*ip->cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  smax = SHRT_MIN;
  smin = SHRT_MAX;

  //  if(strcmp(ip->mediasopUID,ip->sopUID))
  //    printf("Whoa, dude!  I didn't think this could happen!  I'm scared!\n");

  strcpy(uid,ip->sopUID);
  length = strlen(uid);
  for(i=7; i>0; i--) {  // randomize last seven digits
    uid[length-i] += (char)(ran1(&ldum)*10.0);
    if(uid[length-i] > '9')
      uid[length-i] -= 10;
  }
  
  if(!seriesID[seriesno - firstseries]) { // first image for this series
    strcpy(suid,ip->seriesUID);
    length = strlen(suid);
    for(oobs = strlen(suid)-1;(suid[oobs] != '.') && (oobs > -1); oobs--);
    for(oobs--;(suid[oobs] != '.') && (oobs > -1); oobs--);  //second
    for(oobs--;(suid[oobs] != '.') && (oobs > -1); oobs--);  //third
    if(oobs < 7) {
       fprintf(stderr,"error parsing series UID %s, bailing!\n",ip->seriesUID);
      exit(1);
    }
    for(i=7; i>0; i--) {  // randomize last seven digits of appropriate field
      suid[oobs-i] += (char)(ran1(&ldum)*10.0);
      if(suid[oobs-i] > '9')
	suid[oobs-i] -= 10;
    }
    if(strlen(suid)%2) {  // needs a padding space at end to be even
      suid[strlen(suid)+1] = '\0';
      suid[strlen(suid)] = ' ';
    }
    if((seriesID[seriesno - firstseries] = malloc(strlen(suid)+1)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    strcpy(seriesID[seriesno - firstseries], suid);
    //    printf("made new series idea from %s\n",ip->seriesUID);
  }

  for(i=0; i<ip->rows*ip->cols; i++) {
    sbuf[i] = (short)dbuf[i];
    if(sbuf[i] > smax)
      smax = sbuf[i];
    if(sbuf[i] < smin)
      smin = sbuf[i];
  }
  
  done = 0;
  didmyfields = 0;
  if((ifp = fopen(ip->filename,"r")) == NULL) {
    fprintf(stderr,"could not open file %s: ",ip->filename);
    perror("");
    exit(1);
  }
  //  sprintf(filename,"%s%.5s\0",fname,ip->sliceloc);
  if(basename)
    sprintf(filename,"%s%s%03d\0",basename,fname,imageno);
  else 
    sprintf(filename,"%s%03d\0",fname,imageno);
  //  printf("creating file %s based on file %s\n",filename,ip->filename);
  if(destdir[0]) { // might need to make directory
    if(stat(destdir,&statbuf)) {
      printf("output directory %s does not exist, bailing...\n",destdir);
      exit(1);
    }
    sprintf(filename,"%s/%s\0",destdir,ip->name);
    if(stat(filename,&statbuf)) {
      printf("making directory %s\n",filename);
      if(mkdir(filename,0777)) {
	fprintf(stderr,"mkdir failed for directory %s: ",filename);
	perror("");
	exit(1);
      }
    }
    else if(!S_ISDIR(statbuf.st_mode)) {
      fprintf(stderr,"error: %s exists and is not a directory: ",filename);
      perror("");
      exit(1);
    }
    
    sprintf(filename,"%s/%s/%s\0",destdir,ip->name,ip->date);
    if(stat(filename,&statbuf)) {
      printf("making directory %s\n",filename);
      if(mkdir(filename,0777)) {
	fprintf(stderr,"mkdir failed for directory %s: ",filename);
	perror("");
	exit(1);
      }
    }
    else if(!S_ISDIR(statbuf.st_mode)) {
      fprintf(stderr,"error: %s exists and is not a directory: ",filename);
      perror("");
      exit(1);
    }
    if(!resultsdir[0])
      strcpy(resultsdir,filename);
    if(basename)
      sprintf(filename,"%s/%s/%s/%s%s%03d\0",destdir,ip->name,ip->date,basename,fname,imageno);
    else
      sprintf(filename,"%s/%s/%s/%s%03d\0",destdir,ip->name,ip->date,fname,imageno);

  }
  if(!clobber)
    if(!stat(filename,&statbuf)) {
      printf("file %s exists, skipping...\n",filename);
      fclose(ifp);
      free(sbuf);
      return;
    }
  if((ofp = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"could not open %s for DICOM output: ",filename);
    perror("");
  }
  fwrite(cp,1,128,ofp);
  if((cp = malloc(128)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  cplength = 128;
  fwrite("DICM",1,4,ofp);
  fseek(ifp,132L,SEEK_SET);
  for(done=0;!done;) {
    if(fread(&tag,sizeof(struct tag),1,ifp) != 1) {
      perror("read error 1");
      exit(1);
    }
    if(!strncmp(tag.VR,"OB",2) || !strncmp(tag.VR,"OW",2) || !strncmp(tag.VR,"OF",2)
       || !strncmp(tag.VR,"SQ",2) || !strncmp(tag.VR,"UT",2) || !strncmp(tag.VR,"UN",2)) {
      if(fread(&sdum,2,1,ifp) != 1) {
	perror("read error 2");
	exit(1);
      }
      if(sdum) {
	fprintf(stderr,"inconsistency in VR field, should read 0, instead is %0.4X\n",sdum);
	exit(1);
      }
      if(fread(&length,4,1,ifp) != 1) {
	perror("read error 3");
	exit(1);
      }
    }
    else {
      if(fread(&slength,2,1,ifp) != 1) {
	perror("read error 4");
	exit(1);
      }
      length = (int32_t)slength;
    }
    if(!didmyfields) {
      if(tag.group == 0x21)
	if(tag.element == 0x10)
	  {
	    fprintf(stderr,"hey, someone is already using private tag group 0x21, element 0x10!  Bailing!\n");
	    exit(1);
	  }
      if(tag.group > 0x20) {
	if(debug)
	  printf("inserting my dicom fields...\n");
	pofftag.group = 0x21;
	pofftag.element = 0x10;
	pofftag.VR[0] = 'L';
	pofftag.VR[1] = 'O';
	sprintf(poffstr,"Glen Morrell's Kludge Dixon Software\0");
	if(strlen(poffstr) % 2) {
	  poffstr[strlen(poffstr)+1] = '\0';
	  poffstr[strlen(poffstr)] = ' ';
	}
	if(fwrite(&pofftag,sizeof(struct tag),1,ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	pofflen = strlen(poffstr);
	if(fwrite(&pofflen,sizeof(short),1,ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fwrite(poffstr,1,strlen(poffstr),ofp) != pofflen) {
	  perror("write error");
	  exit(1);
	}
	
	// put command line in group 0x21, element 0x1000
	pofftag.element = 0x1000;
	pofftag.VR[0] = 'L';
	pofftag.VR[1] = 'T';
	sprintf(poffstr,"%s\0",cmdline);
	if(strlen(poffstr) % 2) {
	  poffstr[strlen(poffstr)+1] = '\0';
	  poffstr[strlen(poffstr)] = ' ';
	}
	pofflen = strlen(poffstr);
	if(fwrite(&pofftag,sizeof(struct tag),1,ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fwrite(&pofflen,sizeof(short),1,ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fwrite(poffstr,1,strlen(poffstr),ofp) != pofflen) {
	  perror("write error");
	  exit(1);
	}
	
	pofftag.element = 0x1001;  // timestamp in 0x1001
	dtime = time((time_t *)0);
	sprintf(poffstr,"%s",ctime(&dtime));
	poffstr[strlen(poffstr)-1] = '\0';  // get rid of ending \n
	if(strlen(poffstr) % 2) {
	  poffstr[strlen(poffstr)+1] = '\0';
	  poffstr[strlen(poffstr)] = ' ';
	}
	pofflen = strlen(poffstr);
	if(fwrite(&pofftag,sizeof(struct tag),1,ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fwrite(&pofflen,sizeof(short),1,ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fwrite(poffstr,1,strlen(poffstr),ofp) != pofflen) {
	  perror("write error");
	  exit(1);
	}
	
	// put phase offset in group 0x21, element 0x1002
	if(poff) {
	  pofftag.element = 0x1002;
	  sprintf(poffstr,"poff = %g\n",poff);
	  if(strlen(poffstr) % 2) {
	    poffstr[strlen(poffstr)+1] = '\0';
	    poffstr[strlen(poffstr)] = ' ';
	  }
	  pofflen = strlen(poffstr);
	  if(fwrite(&pofftag,sizeof(struct tag),1,ofp) != 1) {
	    perror("write error");
	    exit(1);
	  }
	  if(fwrite(&pofflen,sizeof(short),1,ofp) != 1) {
	    perror("write error");
	    exit(1);
	  }
	  if(fwrite(poffstr,1,strlen(poffstr),ofp) != pofflen) {
	    perror("write error");
	    exit(1);
	  }
	}
	didmyfields++;
      }
    }
    // Deal with SQ lists
    if(!strncmp(tag.VR,"SQ",2)) {
      if(fwrite(&tag,sizeof(struct tag),1,ofp) != 1) {
	perror("write error");
	exit(1);
      }
      sdum = 0;
      if(fwrite(&sdum,sizeof(short),1,ofp) != 1) {
	perror("write error");
	exit(1);
      }
      if(fwrite(&length,sizeof(int32_t),1,ofp) != 1) {
	perror("write error");
	exit(1);
      }
      if(debug) {
	printf("writing tag group %0.4X, element %0.4X, VR %c%c length %ld\n",tag.group,tag.element,tag.VR[0],tag.VR[1],length);
	printf("calling sqloop with length %ld\n",length);
      }
      sqloop(&ifp,&ofp,length);
      continue;
    }
    if(length < 0) {
      printf("whoa, negative length outside of SQ list!  Bailing!\n");
      exit(1);
    }
    if(cplength < (length+1)) {
      if((cp = realloc(cp,length+1)) == NULL) {
	perror("realloc failed");
	exit(1);
      }
    }
    if(fread(cp,1,length,ifp) != length) {
      perror("read error 5");
      exit(1);
    }
    switch(tag.group) {
    case 0x2:
      switch(tag.element) {
      case 0x3:
	strcpy(cp,uid);
	length = strlen(cp);
	break;
      default:
	break;
      }
      break;
    case 0x8:
      switch(tag.element) {
      case 0x18:
	strcpy(cp,uid);
	length = strlen(cp);
	break;
      case 0x103E:  // series description
	strcpy(cp,protstring);
	length = strlen(protstring);
	if(length % 2) {
	  cp[length] = ' ';
	  cp[++length] = '\0';
	}
	break;
      default:
	break;
      }
      break;
    case 0x18:
      switch(tag.element) {
      case 0x1030:
	strcpy(cp,protstring);
	length = strlen(protstring);
	if(length % 2) {
	  cp[length] = ' ';
	  cp[++length] = '\0';
	}
	break;
      default:
	break;
      }
      break;
    case 0x20:
      switch(tag.element) {
      case 0xe:
	strcpy(cp,seriesID[seriesno - firstseries]);
	//	printf("writing series id %s\n",seriesID[seriesno - firstseries]);
	length = strlen(cp);
	break;
      case 0x11: // series number
	sprintf(cp,"%d\0",seriesno);
	length = strlen(cp);
	if(length % 2) {  // must pad out to even number of characters.
	  cp[length] = ' ';
	  cp[length+1] = '\0';
	  length++;
	}
	break;
      case 0x13: // image number
	sprintf(cp,"%d\0",imageno);
	length = strlen(cp);
	if(length % 2) {  // must pad out to even number of characters.
	  cp[length] = ' ';
	  cp[length+1] = '\0';
	  length++;
	}
	break;
      default:
	break;
      }
      break;
    case 0x28:
      switch(tag.element) {
      case 0x10:
	schar.s = ip->rows;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];  // used strncpy here and didn't work because stops at null char
	break;
      case 0x11:
	schar.s = ip->cols;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	break;
      case 0x100:  // bits allocated
	schar.s = sizeof(short) * 8;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	break;
      case 0x101:  // bits stored
	schar.s = sizeof(short) * 8;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	break;
      case 0x102:  // high bit
	schar.s = sizeof(short) * 8 -1;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	break;
      case 0x103:  // pixel representation, seems that 1 = 16 bit signed
	schar.s = 1;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	break;
      case 0x106:  // min pixel value
	schar.s = smin;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	tag.VR[0] = 'S';
	tag.VR[1] = 'S';  // signed short
	break;
      case 0x107:  // max pixel value
	schar.s = smax;
	for(i=0; i<length; i++)
	  cp[i] = schar.c[i];
	tag.VR[0] = 'S';
	tag.VR[1] = 'S';  // signed short
	break;
      case 0x1050:   // window center
	sprintf(cp,"%d\0",(int)((double)(smin + smax)/2.0));
	length = strlen(cp);
	if(length % 2) {  // must pad out to even number of characters.
	  cp[length] = ' ';
	  cp[length+1] = '\0';
	  length++;
	}
	break;
      case 0x1051:   // window width
	sprintf(cp,"%d\0",(int)(smax - smin));
	length = strlen(cp);
	if(length % 2) {  // must pad out to even number of characters.
	  cp[length] = ' ';
	  cp[length+1] = '\0';
	  length++;
	}
	break;
      case 0x1052:
      case 0x1053:
      case 0x1054: // bad fields here, get set in the phase images
	continue;
      default:
	break;
      }
      break;
    case 0x7fe0:
      switch(tag.element) {
      case 0x10:
	length = ip->rows * ip->cols * sizeof(short);
	done++;
	break;
      default:
	break;
      }
      break;
    default:
      break;
    }
    if(debug)
      printf("writing tag group %0.4X, element %0.4X, VR %c%c length %ld\n",tag.group,tag.element,tag.VR[0],tag.VR[1],length);
    if(fwrite(&tag,sizeof(struct tag),1,ofp) != 1) {
      perror("write error");
      exit(1);
    }
    if(!strncmp(tag.VR,"OB",2) || !strncmp(tag.VR,"OW",2) || !strncmp(tag.VR,"OF",2)
       || !strncmp(tag.VR,"SQ",2) || !strncmp(tag.VR,"UT",2) || !strncmp(tag.VR,"UN",2)) {
      sdum = 0;
      if(fwrite(&sdum,sizeof(short),1,ofp) != 1) {
	perror("write error");
	exit(1);
      }
      if(fwrite(&length,sizeof(int32_t),1,ofp) != 1) {
	perror("write error");
	exit(1);
      }
    }
    else {
      sdum = (short)length;
      if(fwrite(&sdum,sizeof(short),1,ofp) != 1) {
	perror("write error");
	exit(1);
      }
    }
    if(!done)
      if(fwrite(cp,1,length,ofp) != length) {
	perror("write error");
	exit(1);
      }
  }
  // Assuming image data is always last item in the DICOM file.
  if(fwrite(sbuf,sizeof(short),ip->rows*ip->cols,ofp) != (ip->rows*ip->cols)) {
    perror("write error");
    exit(1);
  }
  fclose(ifp);
  fclose(ofp);
  free(sbuf);
}

writeRaw(fname,ip,dbuf)
     char *fname;
     struct image *ip;
     double *dbuf;
{
  FILE *fp;
  short *sbuf;
  struct stat statbuf;
  int i;
  char filename[1024];

  if((sbuf = (short *)malloc(sizeof(short)*ip->rows*ip->cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if(destdir[0]) { // might need to make directory
    if(stat(destdir,&statbuf)) {
      printf("output directory %s does not exist, bailing...\n",destdir);
      exit(1);
    }
    sprintf(filename,"%s/%s\0",destdir,ip->name);
    if(stat(filename,&statbuf)) {
      printf("making directory %s\n",filename);
      if(mkdir(filename,0777)) {
	fprintf(stderr,"mkdir failed for directory %s: ",filename);
	perror("");
	exit(1);
      }
    }
    else if(!S_ISDIR(statbuf.st_mode)) {
      fprintf(stderr,"error: %s exists and is not a directory: ",filename);
      perror("");
      exit(1);
    }
    
    sprintf(filename,"%s/%s/%s\0",destdir,ip->name,ip->date);
    if(stat(filename,&statbuf)) {
      printf("making directory %s\n",filename);
      if(mkdir(filename,0777)) {
	fprintf(stderr,"mkdir failed for directory %s: ",filename);
	perror("");
	exit(1);
      }
    }
    else if(!S_ISDIR(statbuf.st_mode)) {
      fprintf(stderr,"error: %s exists and is not a directory: ",filename);
      perror("");
      exit(1);
    }
    if(basename)
      sprintf(filename,"%s/%s/%s/%s%s\0",destdir,ip->name,ip->date,basename,fname);
    else
      sprintf(filename,"%s/%s/%s/%s\0",destdir,ip->name,ip->date,fname);
  }
  else
    if(basename)
      sprintf(filename,"%s%s",basename,fname);
    else
      sprintf(filename,"%s",fname);
  if(!clobber)
    if(!stat(filename,&statbuf)) {
      printf("file %s exists, skipping...\n",filename);
      free(sbuf);
      return;
    }
  if((fp = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"could not open %s: ",filename);
    perror("");
    exit(1);
  }
  for(i=0; i<ip->rows*ip->cols; i++)
    sbuf[i] = (short)dbuf[i];
  /*  fwrite(dbuf,sizeof(double),ip->rows*ip->cols,debug);
  fwrite(sbuf,sizeof(short),ip->rows*ip->cols,debug);
  printf("stopped in rawWrite before writing %s\n",filename);
  exit(0);
  */

  if(fwrite(sbuf,sizeof(short),ip->rows*ip->cols,fp) != ip->rows*ip->cols) {
    perror("write error");
    exit(1);
  }
  fclose(fp);
  free(sbuf);
}

sqloop(ifp,ofp,listlength) 
     FILE **ifp,**ofp;
     int32_t listlength;
{
  struct tag tag;
  int32_t cplength,length;
  unsigned short sdum,slength;
  short group,element;
  char *cp;

  cplength = 0;
  cp = (char *)0;
  if(listlength != -1) {  // easy, SQ list with explicit length
    if((cp = malloc(listlength)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    if(fread(cp,1,listlength,*ifp) != listlength) {
      perror("read error 6");
      exit(1);
    }
    if(fwrite(cp,1,listlength,*ofp) != listlength) {
      perror("write error");
      exit(1);
    }
    if(debug)
      printf("in sqloop, explicit length, wrote %ld bytes.\n",listlength);
    free(cp);
    return;
  }
  // listlength == -1, delimited list
  for(;;) {
  kludge:
    if(fread(&sdum,2,1,*ifp) != 1) {
      perror("read error 7");
      exit(1);
    }
    if(debug) 
      printf("item: bytes %0.4X, ",sdum);
    if(sdum != 0xFFFE) {  // each item is supposed to start with 0xFFFE
      fprintf(stderr,"in writeDicom, item in SQ list did not start with 0xFFFE.  Bailing!\n");
      exit(1);
    }
    if(fwrite(&sdum,2,1,*ofp) != 1) {
      perror("write error");
      exit(1);
    }
    if(fread(&sdum,2,1,*ifp) != 1) {
      perror("read error 8");
      exit(1);
    }
    if(debug)
      printf("%0.4X length ",sdum);
    if(fwrite(&sdum,2,1,*ofp) != 1) {
      perror("write error");
      exit(1);
    }
    if(fread(&length,4,1,*ifp) != 1) {
      perror("read error 9");
      exit(1);
    }
    if(fwrite(&length,4,1,*ofp) != 1) {
      perror("write error");
      exit(1);
    }
    if(debug)
      printf("%ld\n",length);
    if(sdum == 0xE0DD) { // code for end of list, followed by zero length (4 bytes)
      if(length) {
	fprintf(stderr,"unexpected non-zero length %ld after sequence delimiter, bailing...\n",length);
	exit(1);
      }
      //      if(fwrite(&length,4,1,*ofp) != 1) {
      //	perror("write error");
      //	exit(1);
      //      }
      if(debug)
	printf("end of list, returning from sqloop\n");
      return;
    }
    if(sdum != 0xE000) { // code for beginning of new item
      fprintf(stderr,"unexpected group/element 0xFFFE, 0x%.X\n",sdum);
      exit(1);
    }
    if(length != -1) {  // item has explicit length
      if(length > cplength) {
	if(cp) {
	  if((cp = realloc(cp,length)) == NULL) {
	    perror("realloc failed");
	    exit(1);
	  }
	} else {
	  if((cp = malloc(length)) == NULL) {
	    perror("malloc failed");
	    exit(1);
	  }
	}
      }
      if(fread(cp,1,length,*ifp) != length) {
	perror("read error 10");
	exit(1);
      }
      if(fwrite(cp,1,length,*ofp) != length) {
	perror("write error");
	exit(1);
      }
      if(debug)
	printf("wrote item with explicit length %ld\n",length);
      continue;
    }
    // item with non-specified length

    for(;;){  // loop through data elements in this item looking for nested SQ lists
      if(fread(&tag.group,2,1,*ifp) != 1) {
	perror("read error 11");
	exit(1);
      }
      if(fread(&tag.element,2,1,*ifp) != 1) {
	perror("read error 12");
	exit(1);
      }
      if((tag.group == 0xFFFE) && (tag.element == 0xE00D)) { // end of item
	if(fwrite(&tag.group,2,1,*ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fwrite(&tag.element,2,1,*ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fread(&length,4,1,*ifp) != 1) {
	  perror("read error 13");
	  exit(1);
	}
	if(length) {
	  fprintf(stderr,"got non-zero length for end of item, bailing...\n");
	  exit(1);
	}
	if(fwrite(&length,4,1,*ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(debug)
	  printf("in sqlist, group = %0.4X element = %0.4X, end of item\n",tag.group,tag.element);
	goto kludge;
      }
      if(fread(&tag.VR,2,1,*ifp) != 1) {
	perror("read error 14");
	exit(1);
      }
      if(fwrite(&tag,sizeof(struct tag),1,*ofp) != 1) {
	perror("write error");
	exit(1);
      }
      if(!strncmp(tag.VR,"OB",2) || !strncmp(tag.VR,"OW",2) || !strncmp(tag.VR,"OF",2)
	 || !strncmp(tag.VR,"SQ",2) || !strncmp(tag.VR,"UT",2) || !strncmp(tag.VR,"UN",2)){
	
	if(fread(&sdum,2,1,*ifp) != 1) {
	  perror("read error 15");
	  exit(1);
	}
	if(sdum) {
	  fprintf(stderr,"inconsistency in VR field, should read 0, instead is %0.4X\n",sdum);
	  exit(1);
	}
	if(fwrite(&sdum,2,1,*ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	if(fread(&length,4,1,*ifp) != 1) {
	  perror("read error 16");
	  exit(1);
	}
	if(fwrite(&length,4,1,*ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
      }
      else {
	if(fread(&slength,2,1,*ifp) != 1) {
	  perror("read error 17");
	  exit(1);
	}
	if(fwrite(&slength,2,1,*ofp) != 1) {
	  perror("write error");
	  exit(1);
	}
	length = (int32_t)slength;
      }
      if(debug)
	printf("writing tag group %0.4X, element %0.4X, VR %c%c length %ld\n",tag.group,tag.element,tag.VR[0],tag.VR[1],length);
      if(!strncmp(tag.VR,"SQ",2)) {
	if(debug)
	  printf("calling sqloop with length %ld\n",length);
	sqloop(ifp,ofp,length);
	continue;
      }
      if(length > cplength) {
	if(cp) {
	  if((cp = realloc(cp,length)) == NULL) {
	    perror("realloc failed");
	    exit(1);
	  }
	} else {
	  if((cp = malloc(length)) == NULL) {
	    perror("malloc failed");
	    exit(1);
	  }
	}
      }
      if(fread(cp,1,length,*ifp) != length) {
	perror("read error 18");
	exit(1);
      }
      if(fwrite(cp,1,length,*ofp) != length) {
	perror("write error");
	exit(1);
      }
    }
  }
  free(cp);
}
