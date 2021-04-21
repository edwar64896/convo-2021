/*

procedure for creating correct header files:

1. Inspect the WAV file to be used as an impulse response:
	It should have 4 channels, be PCM and have 24 bit resolution:
2. Convert to headerless raw using sndfile-convert
3. convert to a .h file using xxd -i
4. include this file in the header set

*/

/*
typedef wf_count_t  (*sf_vio_get_filelen) (void *user_data) ;
typedef sf_count_t  (*sf_vio_seek)        (sf_count_t offset, int whence, void *user_data) ;
typedef sf_count_t  (*sf_vio_read)        (void *ptr, sf_count_t count, void *user_data) ;
typedef sf_count_t  (*sf_vio_write)       (const void *ptr, sf_count_t count, void *user_data) ;
typedef sf_count_t  (*sf_vio_tell)        (void *user_data) ;
*/


/* these files are to be used as callbacks for the vio functions in libsndfile */

wf_count_t fir1_get_filelen(void *user_data) {
}

wf_count_t fir2_get_filelen(void *user_data) {
}

sf_count_t fir1_vio_seek(sf_count_t offset, int whence, void *user_data) {
}

sf_count_t fir2_vio_seek(sf_count_t offset, int whence, void *user_data) {
}

sf_count_t fir1_vio_read(void *ptr, sf_count_t count, void *user_data) {
}

sf_count_t fir2_vio_read(void *ptr, sf_count_t count, void *user_data) {
}

sf_count_t fir1_vio_write (const void *ptr, sf_count_t count, void *user_data) {
}

sf_count_t fir2_vio_write (const void *ptr, sf_count_t count, void *user_data) {
}

sf_count_t fir1_vio_tell (void *user_data){
}

sf_count_t fir2_vio_tell (void *user_data){
}
