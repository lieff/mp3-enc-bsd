
/*
  Implements the pseudocode of page 98 of the IS

  1: 1027337 => 1.819941
  2: 1031950 => 1.828113
  3: 1027337 => 1.819941
  4: 3052971 => 5.408377
  5: 1031950 => 1.828113
  6: 3052510 => 5.407560
  
  7: 23114009 => 40.946761
  8: 23110868 => 40.941196

  */
static inline int
HuffmanCode( int table_select, int x, int y, BF_PartHolder **pph )
{
    unsigned int signx = 0, signy = 0, linbitsx, linbitsy, linbits, xlen, ylen, idx;
    unsigned int code, ext;
    int          cbits, xbits;
    struct huffcodetab *h;

    cbits = 0;
    xbits = 0;
    code  = 0;
    ext   = 0;
    
    if ( table_select == 0 )
	return 0;
    
#if 0  /* 4% cpu wasted here! */
    signx = abs_and_sign( &x );
    signy = abs_and_sign( &y );
#else
    if ( x < 0 ) {
      x     = -x;
      signx = 1;
    }
    if ( y < 0 ) {
      y     = -y;
      signy = 1;
    }
#endif
    h = &(ht[table_select]);
    xlen = h->xlen;
    ylen = h->ylen;
    linbits = h->linbits;
    linbitsx = linbitsy = 0;

    if ( table_select > 15 )
    { /* ESC-table is used */
	if ( x > 14 )
	{
	    linbitsx = x - 15;
//	    assert( linbitsx <= h->linmax );
	    x = 15;
	}
	if ( y > 14 )
	{
	    linbitsy = y - 15;
//	    assert( linbitsy <= h->linmax );
	    y = 15;
	}
	idx = (x * ylen) + y;
	code = h->table[idx];
	cbits = h->hlen[ idx ];
	if ( x > 14 )
	{
	    ext |= linbitsx;
	    xbits += linbits;
	}
	if ( x != 0 )
	{
	    ext <<= 1;
	    ext |= signx;
	    xbits += 1;
	}
	if ( y > 14 )
	{
	    ext <<= linbits;
	    ext |= linbitsy;
	    xbits += linbits;
	}
	if ( y != 0 )
	{
	    ext <<= 1;
	    ext |= signy;
	    xbits += 1;
	}
    }
    else
    { /* No ESC-words */
	idx = (x * ylen) + y;
	code = h->table[idx];
	cbits += h->hlen[ idx ];
	if ( x != 0 )
	{
	    code <<= 1;
	    code |= signx;
	    cbits += 1;
	}
	if ( y != 0 )
	{
	    code <<= 1;
	    code |= signy;
	    cbits += 1;
	}
    }
//    assert( cbits <= 32 );
//    assert( xbits <= 32 );

    if (pph) {
#if 0
      *pph = BF_addEntry( *pph,  code, cbits );
      *pph = BF_addEntry( *pph,  ext, xbits );
#else

      BF_BitstreamElement myElement;
      BF_PartHolder *thePH = *pph;
      
      myElement.value  = code;
      myElement.length = cbits;
      if ( cbits )
	thePH = BF_addElement( thePH, &myElement );
      
      myElement.value  = ext;
      myElement.length = xbits;
      if ( xbits )
	thePH = BF_addElement( thePH, &myElement );
      
      if ( xbits || xbits )
	*pph = thePH;
      
#endif
    }
    
    return cbits + xbits;
}
