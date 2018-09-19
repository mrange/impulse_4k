//--------------------------------------------------------------------------//
// iq / rgba  .  tiny codes  .  2008/2015                                   //
//--------------------------------------------------------------------------//

#define WIN32_LEAN_AND_MEAN
#define WIN32_EXTRA_LEAN
#include <windows.h>
#include <GL/gl.h>
#include <math.h>
#include "config.h"
#include "ext.h"
#include "shader.h"
#include "fp.h"

//=================================================================================================================

static int   fsid;

int intro_init( void )
{
    if( !EXT_Init() )
        return 0;

    char * vsh = getVertexShader();
    char * fsh = getFragmentShader();  

    int vsid = oglCreateShaderProgramv( GL_VERTEX_SHADER,   1, &vsh );
        fsid = oglCreateShaderProgramv( GL_FRAGMENT_SHADER, 1, &fsh );
 
    unsigned int pid;
    oglGenProgramPipelines(1, &pid);
    oglBindProgramPipeline(pid);
    oglUseProgramStages(pid, GL_VERTEX_SHADER_BIT, vsid);
    oglUseProgramStages(pid, GL_FRAGMENT_SHADER_BIT, fsid);

    #ifdef DEBUG
        int     result;
        char    info[1536];
        oglGetProgramiv( vsid, GL_LINK_STATUS, &result ); oglGetProgramInfoLog( vsid, 1024, NULL, (char *)info );
        // if( !result ) DebugBreak();
        oglGetProgramiv( fsid, GL_LINK_STATUS, &result ); oglGetProgramInfoLog( fsid, 1024, NULL, (char *)info ); 
        // if( !result ) DebugBreak();
        oglGetProgramiv( pid,  GL_LINK_STATUS, &result ); oglGetProgramInfoLog( pid,  1024, NULL, (char *)info ); 
        // if( !result ) DebugBreak();
    #endif

    return 1;
}

//=================================================================================================================

static float fparams[4*1];

void intro_do( long time, int xres, int yres )
{
    //--- update parameters -----------------------------------------

    const float t  = 0.001f*(float)time;

    // camera position
    fparams[ 0] = t;
    fparams[ 1] = xres;
    fparams[ 2] = yres;
    fparams[ 3] = 0.0f;

    //--- render -----------------------------------------

    oglProgramUniform4fv( fsid, 0, 1, fparams );

    glRects( -1, -1, 1, 1 );
}