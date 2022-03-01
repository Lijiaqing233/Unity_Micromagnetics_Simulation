using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public struct double3
{
    public double x, y, z;
}

public struct float3
{
    public float x, y, z;
}




public struct Particle
{//8*9 = 72
    public double qxx, qyy, qzz, qxy, qxz, qyz;
    public double li, lj, lk;
    //24*8 = 192
    public double3 vec, b, v, K1, K2, K3, K4, position;//264
    
}

public struct Particlef
{//8*9 = 72
    public float qxx, qyy, qzz, qxy, qxz, qyz;
    public float li, lj, lk;
    //24*8 = 192
    public float3 vec, b, v, K1, K2, K3, K4, position;//264

}


public struct Parameter
{ //Number of conductors
    public int nx, ny, nz;
    public int n;//nx * ny * nz;
    public double dx, dy, dz;
    //The magnitude of the external magnetic field
    public double hxc, hyc, hzc, tmax, ipr;
    public int ic;
    public double al, dt, ft;
    public double smx, smy, smz;
    /**********-----Internal fixed parameters------**************/
    public double pi;
    // parameters
    public double ku, m, aa;
    // parameter for dynamic cal
    public double gamma;
}


public struct Parameterf
{ //Number of conductors
    public int nx, ny, nz;
    public int n;//nx * ny * nz;
    public float dx, dy, dz;
    //The magnitude of the external magnetic field
    public float hxc, hyc, hzc, tmax, ipr;
    public int ic;
    public float al, dt, ft;
    public float smx, smy, smz;
    /**********-----Internal fixed parameters------**************/
    public float pi;
    // parameters
    public float ku, m, aa;
    // parameter for dynamic cal
    public float gamma;
}



public class SimBase : MonoBehaviour
{
    /******************************************************/
    public int step = 0;
    public Mesh myvector;

   
    //private Particle tempParticle = new Particle();

    public void init(Parameter[] _Parameter, Particle[] particleList,int nx, int ny , int nz , double hxc, double hyc,
        double hzc, double dx, double dy, double dz, double dt, double alpha,bool ic)
    {
        _Parameter[0].nx = nx;     //导体的数目
        _Parameter[0].ny = ny;
        _Parameter[0].nz = nz;

        _Parameter[0].n = nx * ny * nz;//nx * ny * nz;

        _Parameter[0].dx = (double)(dx * 1e-7);
        _Parameter[0].dy = (double)(dy * 1e-7);
        _Parameter[0].dz = (double)(dz * 1e-7);

        _Parameter[0].hxc = hxc;  //外部磁场的大小
        _Parameter[0].hyc = hyc;
        _Parameter[0].hzc = hzc;
        _Parameter[0].tmax = (double)3e-9;
        _Parameter[0].ipr = 40;
        _Parameter[0].ic = 0;
        _Parameter[0].al = alpha;
        _Parameter[0].dt = (double)(dt * 1e-9);
        _Parameter[0].ft = 0;

        _Parameter[0].smx = 0;
        _Parameter[0].smy = 0;
        _Parameter[0].smz = 0;

        //内部参数
        _Parameter[0].pi = (double)3.141592653589793238462643;
        // parameters
        _Parameter[0].ku = (double)1.85e4;
        _Parameter[0].m = (double)370e0;
        _Parameter[0].aa = (double)0.5e-6;  //0.5e-6正确
                                    // parameter for dynamic cal
        _Parameter[0].gamma = (double)1.76e7; //1.76e7;
        double hk = (double)(2e0 * _Parameter[0].ku / _Parameter[0].m);  //自旋磁场的参数 

        //particleList.Add(tempParticle);
        int temp_index = 0;

        if (ic == true)
        {
            for (int k = 1; k <= _Parameter[0].nz; k++)
            {
                for (int j = 1; j <= _Parameter[0].ny; j++)
                {
                    for (int i = 1; i <= _Parameter[0].nx; i++)
                    {
                        temp_index++;
                        
                        particleList[temp_index].position.x = i * 5 - 10;
                        particleList[temp_index].position.y = j * 5 - 10;
                        particleList[temp_index].position.z = k * 5 - 10;

                        particleList[temp_index].vec.x = 0;
                        particleList[temp_index].vec.y = 0;
                        particleList[temp_index].vec.z = 1;

                        particleList[temp_index].qxx = 0;
                        particleList[temp_index].qyy = 0;
                        particleList[temp_index].qzz = 0;
                        particleList[temp_index].qxy = 0;
                        particleList[temp_index].qxz = 0;
                        particleList[temp_index].qyz = 0;

                        //----------FFT静磁场参数的初始化
                        particleList[temp_index].li = i;
                        particleList[temp_index].lj = j;
                        particleList[temp_index].lk = k;

                        //Debug.Log(k);
                        //particleList.Add(tempParticle);

                    }
                }
            }
        }

        else
        {
            for (int i = 1; i <= _Parameter[0].n; i++)
            {
                 particleList[i].qxx = 0;
                 particleList[i].qyy = 0;
                 particleList[i].qzz = 0;
                 particleList[i].qxy = 0;
                 particleList[i].qxz = 0;
                 particleList[i].qyz = 0;

            }
        }

        //----------FFT静磁场参数的初始化
        
        double s, x, y, yy, z, zz;
        int kc, jc;
        double r;
        
        for (int kk = 0; kk <= 1; kk++)
        {
            for (int jj = 0; jj <= 1; jj++)
            {
                for (int ii = 0; ii <= 1; ii++)
                {
                    s = (double)(_Parameter[0].m * System.Math.Pow(-1e0, (ii + jj + kk)));
                    for (int k = 1; k <= nz ; k++)
                    {
                        z = (double)((k + kk - 1.5e0) * _Parameter[0].dz);
                        zz = z * z;
                        kc = (k - 1) * _Parameter[0].ny;
                        for (int j = 1; j <= _Parameter[0].ny ; j++)
                        {
                            y = (double)((j + jj - 1.5e0) * _Parameter[0].dy);
                            yy = y * y;
                            jc = (kc + j - 1) * _Parameter[0].nx;
                            for (int i = 1; i <= _Parameter[0].nx; i++)
                            {
                                x = (double)((i + ii - 1.5e0) * _Parameter[0].dx);
                                r = (double)(System.Math.Sqrt(x * x + yy + zz));
                                _Parameter[0].ic = jc + i;
                                particleList[_Parameter[0].ic].qxx += (double)(s * System.Math.Atan(y * z / (r * x)));
                                particleList[_Parameter[0].ic].qyy += (double)(s * System.Math.Atan(x * z / (r * y)));
                                particleList[_Parameter[0].ic].qzz += (double)(s * System.Math.Atan(x * y / (r * z)));
                                particleList[_Parameter[0].ic].qxy += (double)(-s * System.Math.Log(System.Math.Abs(z + r)));
                                particleList[_Parameter[0].ic].qxz += (double)(-s * System.Math.Log(System.Math.Abs(y + r)));
                                particleList[_Parameter[0].ic].qyz += (double)(-s * System.Math.Log(System.Math.Abs(x + r)));

                                //Debug.Log("y：" + x + y + z);
                            } 
                        }
                    }
                }
            }
        }

    }


    public void initFloat(Parameterf[] _Parameter, Particlef[] particleList, int nx, int ny, int nz, float hxc, float hyc,
       float hzc, float dx, float dy, float dz, float dt, float alpha, bool ic)
    {
        _Parameter[0].nx = nx;     //导体的数目
        _Parameter[0].ny = ny;
        _Parameter[0].nz = nz;

        _Parameter[0].n = nx * ny * nz;//nx * ny * nz;

        _Parameter[0].dx = (float)(dx * 1e-7);
        _Parameter[0].dy = (float)(dy * 1e-7);
        _Parameter[0].dz = (float)(dz * 1e-7);

        _Parameter[0].hxc = hxc;  //外部磁场的大小
        _Parameter[0].hyc = hyc;
        _Parameter[0].hzc = hzc;
        _Parameter[0].tmax = (float)3e-9;
        _Parameter[0].ipr = 40;
        _Parameter[0].ic = 0;
        _Parameter[0].al = alpha;
        _Parameter[0].dt = (float)(dt * 1e-9);
        _Parameter[0].ft = 0;

        _Parameter[0].smx = 0;
        _Parameter[0].smy = 0;
        _Parameter[0].smz = 0;

        //内部参数
        _Parameter[0].pi = (float)3.141592653589793238462643;
        // parameters
        _Parameter[0].ku = (float)1.85e4;
        _Parameter[0].m = (float)370e0;
        _Parameter[0].aa = (float)0.5e-6;  //0.5e-6正确
                                            // parameter for dynamic cal
        _Parameter[0].gamma = (float)1.76e7; //1.76e7;
        float hk = (float)(2e0 * _Parameter[0].ku / _Parameter[0].m);  //自旋磁场的参数 

        //particleList.Add(tempParticle);
        int temp_index = 0;

        if (ic == true)
        {
            for (int k = 1; k <= _Parameter[0].nz; k++)
            {
                for (int j = 1; j <= _Parameter[0].ny; j++)
                {
                    for (int i = 1; i <= _Parameter[0].nx; i++)
                    {
                        temp_index++;

                        particleList[temp_index].position.x = i * 5 - 10;
                        particleList[temp_index].position.y = j * 5 - 10;
                        particleList[temp_index].position.z = k * 5 - 10;

                        particleList[temp_index].vec.x = 0;
                        particleList[temp_index].vec.y = 0;
                        particleList[temp_index].vec.z = 1;

                        particleList[temp_index].qxx = 0;
                        particleList[temp_index].qyy = 0;
                        particleList[temp_index].qzz = 0;
                        particleList[temp_index].qxy = 0;
                        particleList[temp_index].qxz = 0;
                        particleList[temp_index].qyz = 0;

                        //----------FFT静磁场参数的初始化
                        particleList[temp_index].li = i;
                        particleList[temp_index].lj = j;
                        particleList[temp_index].lk = k;

                        //Debug.Log(k);
                        //particleList.Add(tempParticle);

                    }
                }
            }
        }

        else
        {
            for (int i = 1; i <= _Parameter[0].n; i++)
            {
                particleList[i].qxx = 0;
                particleList[i].qyy = 0;
                particleList[i].qzz = 0;
                particleList[i].qxy = 0;
                particleList[i].qxz = 0;
                particleList[i].qyz = 0;

            }
        }

        //----------FFT静磁场参数的初始化

        float s, x, y, yy, z, zz;
        int kc, jc;
        float r;

        for (int kk = 0; kk <= 1; kk++)
        {
            for (int jj = 0; jj <= 1; jj++)
            {
                for (int ii = 0; ii <= 1; ii++)
                {
                    s = (float)(_Parameter[0].m * System.Math.Pow(-1e0, (ii + jj + kk)));
                    for (int k = 1; k <= nz; k++)
                    {
                        z = (float)((k + kk - 1.5e0) * _Parameter[0].dz);
                        zz = z * z;
                        kc = (k - 1) * _Parameter[0].ny;
                        for (int j = 1; j <= _Parameter[0].ny; j++)
                        {
                            y = (float)((j + jj - 1.5e0) * _Parameter[0].dy);
                            yy = y * y;
                            jc = (kc + j - 1) * _Parameter[0].nx;
                            for (int i = 1; i <= _Parameter[0].nx; i++)
                            {
                                x = (float)((i + ii - 1.5e0) * _Parameter[0].dx);
                                r = (float)(System.Math.Sqrt(x * x + yy + zz));
                                _Parameter[0].ic = jc + i;
                                particleList[_Parameter[0].ic].qxx += (float)(s * System.Math.Atan(y * z / (r * x)));
                                particleList[_Parameter[0].ic].qyy += (float)(s * System.Math.Atan(x * z / (r * y)));
                                particleList[_Parameter[0].ic].qzz += (float)(s * System.Math.Atan(x * y / (r * z)));
                                particleList[_Parameter[0].ic].qxy += (float)(-s * System.Math.Log(System.Math.Abs(z + r)));
                                particleList[_Parameter[0].ic].qxz += (float)(-s * System.Math.Log(System.Math.Abs(y + r)));
                                particleList[_Parameter[0].ic].qyz += (float)(-s * System.Math.Log(System.Math.Abs(x + r)));

                                //Debug.Log("y：" + x + y + z);
                            }
                        }
                    }
                }
            }
        }

    }



    public void LLG (Parameter[] _Parameter, Particle[] particleList, int flag)    //计算llg方程
   {
       int n = _Parameter[0].nx * _Parameter[0].ny * _Parameter[0].nz;
       double3[] hd = new double3[n + 1];
      hd = Demagn(hd, _Parameter, particleList); //静磁场

       // Debug.Log(hd[5].x+"y:"+ hd[5].y + "z:"+hd[5].z);

       //交换磁场计算
       double dx2 = (double)(2e0 * _Parameter[0].aa / _Parameter[0].m / (_Parameter[0].dx * _Parameter[0].dx));
       double dy2 = (double)(2e0 * _Parameter[0].aa / _Parameter[0].m / (_Parameter[0].dy * _Parameter[0].dy));
       double dz2 = (double)(2e0 * _Parameter[0].aa / _Parameter[0].m / (_Parameter[0].dz * _Parameter[0].dz));
       double dyn = (double)(_Parameter[0].gamma * _Parameter[0].dt / (1e0 + _Parameter[0].al * _Parameter[0].al));
       _Parameter[0].ft = (double)0e0;
       int im, ip;
       int jm, jp;
       int km, kp;
       double dmx, dmy, dmz;
       double x, y, z;
       double xxm, xxp, xym, xyp, xzm, xzp;
       double yxm, yxp, yym, yyp, yzm, yzp;
       double zxm, zxp, zym, zyp, zzm, zzp;
       double hxe, hye, hze, hm, tx, ty, tz;
       int m2 = _Parameter[0].nx * _Parameter[0].ny;
       double hk = (double)(2e0 * _Parameter[0].ku / _Parameter[0].m);

       for (int i = 1; i <= n; i++)
       {
           double i1 = particleList[i].li;
           double j1 = particleList[i].lj;
           double k1 = particleList[i].lk;
           //寻找i周围的物体
           im = i - 1;
           ip = i + 1;
           if (i1 <= 1)
               im = i;
           if (i1 >= _Parameter[0].nx)
               ip = i;

           jm = i - _Parameter[0].nx;
           jp = i + _Parameter[0].nx;
           if (j1 <= 1)
               jm = i;
           if (j1 >= _Parameter[0].ny)
           {
               jp = i;
           }

           km = i - m2;
           kp = i + m2;
           if (k1 <= 1)
               km = i;

            if (k1 >= _Parameter[0].nz)
               kp = i;
            
            x = particleList[i].vec.x;
           xxm = particleList[im].vec.x;  //右边
           xxp = particleList[ip].vec.x;  //左边
           xym = particleList[jm].vec.x;   //后面
           xyp = particleList[jp].vec.x;   //前面
           xzm = particleList[km].vec.x;  //下面
            
           xzp = particleList[kp].vec.x;  //上面
           y = particleList[i].vec.y;
           yxm = particleList[im].vec.y;
           yxp = particleList[ip].vec.y;
           yym = particleList[jm].vec.y;
           yyp = particleList[jp].vec.y;
           yzm = particleList[km].vec.y;
           yzp = particleList[kp].vec.y;
           z = particleList[i].vec.z;
           zxm = particleList[im].vec.z;
           zxp = particleList[ip].vec.z;
           zym = particleList[jm].vec.z;
           zyp = particleList[jp].vec.z;
           zzm = particleList[km].vec.z;
           zzp = particleList[kp].vec.z;

           dmx = (double)((xxp - 2e0 * x + xxm) * dx2
                + (xyp - 2e0 * x + xym) * dy2
                + (xzp - 2e0 * x + xzm) * dz2);
           dmy = (double)((yxp - 2e0 * y + yxm) * dx2
                + (yyp - 2e0 * y + yym) * dy2
                + (yzp - 2e0 * y + yzm) * dz2);
           dmz = (double)((zxp - 2e0 * z + zxm) * dx2
                + (zyp - 2e0 * z + zym) * dy2
                + (zzp - 2e0 * z + zzm) * dz2); 

           hxe = dmx + hd[i].x + _Parameter[0].hxc;
           hye = dmy + hd[i].y + _Parameter[0].hyc;
           hze = dmz + hd[i].z + _Parameter[0].hzc + hk * z;
           hm = x * hxe + y * hye + z * hze;
           tx = y * hze - z * hye;
           ty = z * hxe - x * hze;
           tz = x * hye - y * hxe;
           _Parameter[0].ft = _Parameter[0].ft + tx * tx + ty * ty + tz * tz;

          //  Debug.Log("FT1111: " + _Parameter[0].ft);

            if (flag == 1)
            {

                particleList[i].K1.x = (-tx - _Parameter[0].al * (x * hm - hxe)) * dyn;
                particleList[i].K1.y = (-ty - _Parameter[0].al * (y * hm - hye)) * dyn;
                particleList[i].K1.z = (-tz - _Parameter[0].al * (z * hm - hze)) * dyn;

            }
            else if (flag == 2)
            {

                particleList[i].K2.x = (-tx - _Parameter[0].al * (x * hm - hxe)) * dyn;
                particleList[i].K2.y = (-ty - _Parameter[0].al * (y * hm - hye)) * dyn;
                particleList[i].K2.z = (-tz - _Parameter[0].al * (z * hm - hze)) * dyn;

            }
            else if (flag == 3)
            {

                particleList[i].K3.x = (-tx - _Parameter[0].al * (x * hm - hxe)) * dyn;
                particleList[i].K3.y = (-ty - _Parameter[0].al * (y * hm - hye)) * dyn;
                particleList[i].K3.z = (-tz - _Parameter[0].al * (z * hm - hze)) * dyn;

            }
            else if (flag == 4)
            {

                particleList[i].K4.x = (-tx - _Parameter[0].al * (x * hm - hxe)) * dyn;
                particleList[i].K4.y = (-ty - _Parameter[0].al * (y * hm - hye)) * dyn;
                particleList[i].K4.z = (-tz - _Parameter[0].al * (z * hm - hze)) * dyn;

            }
        }
       //交换磁场计算
     
   }


   public double3[] Demagn(double3[] hd, Parameter[] _Parameter, Particle[] particleList)  //计算由于静磁场引起的损耗
   {
       int ip;
       double id, jd, kd;
       double i, j, k;
       double x, y, z;
       double sij = 0;
       double sik = 0;
       double sjk = 0;
       //磁场损耗
       int n = _Parameter[0].nx * _Parameter[0].ny * _Parameter[0].nz;

       for (int tmp = 1; tmp < n + 1; tmp++)
       {
           hd[tmp].x = 0;
           hd[tmp].y = 0;
           hd[tmp].z = 0;
       }
        
       
       for (int ib = 1; ib < n + 1; ib++)
       {
           id = particleList[ib].li;
           jd = particleList[ib].lj;
           kd = particleList[ib].lk;

           x = particleList[ib].vec.x;
           y = particleList[ib].vec.y;
           z = particleList[ib].vec.z;
            
           for (int ii = 1; ii < n + 1; ii++)
           {
               i = particleList[ii].li;
               j = particleList[ii].lj;
               k = particleList[ii].lk;
                double Ks;// = System.Math.Sign(kd - k);
                double Js;// = System.Math.Sign(jd - j);
                double Is;// = System.Math.Sign(id - i);
                if ((kd - k) < 0) { Ks = -1; }
                else { Ks = 1; }
                if ((jd - j) < 0) { Js = -1; }
                else { Js = 1; }
                if ((id - i) < 0) { Is = -1; }
                else { Is = 1; }


                ip = (int)((System.Math.Abs(kd - k) * _Parameter[0].ny + System.Math.Abs(jd - j)) * _Parameter[0].nx + System.Math.Abs(id - i)) + 1;
               sij = Is * Js * particleList[ip].qxy;
               sik = Is * Ks * particleList[ip].qxz;
               sjk = Js * Ks * particleList[ip].qyz;

               hd[ii].x += particleList[ip].qxx * x + sij * y + sik * z;
               hd[ii].y += sij * x + particleList[ip].qyy * y + sjk * z;
               hd[ii].z += sik * x + sjk * y + particleList[ip].qzz * z;
           }
       }
       return hd;
   }

    //ok
    public void vadd(double f, Parameter[] _Parameter, Particle[] particleList, int flag)
    {
        int n = _Parameter[0].n;
        if (flag == 1)
        {
            for (int i = 1; i <= n; i++)
            {
                particleList[i].vec.x = particleList[i].v.x + f * particleList[i].K1.x;
                particleList[i].vec.y = particleList[i].v.y + f * particleList[i].K1.y;
                particleList[i].vec.z = particleList[i].v.z + f * particleList[i].K1.z;


            }
        }
        else if (flag == 2) 
        {
            for (int i = 1; i <= n; i++)
            {
                particleList[i].vec.x = particleList[i].v.x + f * particleList[i].K2.x;
                particleList[i].vec.y = particleList[i].v.y + f * particleList[i].K2.y;
                particleList[i].vec.z = particleList[i].v.z + f * particleList[i].K2.z;


            }
        }
        else if (flag == 3)
        {
            for (int i = 1; i <= n; i++)
            {
                particleList[i].vec.x = particleList[i].v.x + f * particleList[i].K3.x;
                particleList[i].vec.y = particleList[i].v.y + f * particleList[i].K3.y;
                particleList[i].vec.z = particleList[i].v.z + f * particleList[i].K3.z;


            }
        }
        else if (flag == 4)
        {
            for (int i = 1; i <= n; i++)
            {
                particleList[i].vec.x = particleList[i].v.x + f * particleList[i].K4.x;
                particleList[i].vec.y = particleList[i].v.y + f * particleList[i].K4.y;
                particleList[i].vec.z = particleList[i].v.z + f * particleList[i].K4.z;


            }
        }

        for (int i = 1; i <= n; i++)
        {
            double vl = (double)(System.Math.Sqrt(particleList[i].vec.x * particleList[i].vec.x + particleList[i].vec.y * particleList[i].vec.y + particleList[i].vec.z * particleList[i].vec.z));
            if (vl == 0)
                vl = 1;

            particleList[i].vec.x = particleList[i].vec.x / vl;
            particleList[i].vec.y = particleList[i].vec.y / vl;
            particleList[i].vec.z = particleList[i].vec.z / vl;

        }


    }


   public void vadd4(Parameter[] _Parameter, Particle[] particleList)  //ok
   {
       int n = _Parameter[0].nx * _Parameter[0].ny * _Parameter[0].nz;
       for (int i = 1; i < n + 1; i++)
       {
           particleList[i].K4.x = (double)((double)(particleList[i].K1.x + 2e0 * particleList[i].K2.x + 2e0 * particleList[i].K3.x + particleList[i].K4.x) / 6e0);
           particleList[i].K4.y = (double)((double)(particleList[i].K1.y + 2e0 * particleList[i].K2.y + 2e0 * particleList[i].K3.y + particleList[i].K4.y) / 6e0);
           particleList[i].K4.z = (double)((double)(particleList[i].K1.z + 2e0 * particleList[i].K2.z + 2e0 * particleList[i].K3.z + particleList[i].K4.z) / 6e0);
       }

   }

    void vectorCopy(Parameter[] _Parameter, Particle[] particleList) {
    for (int i = 1; i <= _Parameter[0].n;i++) {
        particleList[i].v = particleList[i].vec;
       // cout << "qxx:" << particleList[i].v.x << "     qyy: " << particleList[i].v.y << "    qzz:" << particleList[i].v.z << endl;

    }       
}

   public void rk4(Parameter[] _Parameter, Particle[] particleList)  //4阶龙格库塔方法解算器
   {
        _Parameter[0].ft = 0;
        int n = _Parameter[0].n;
        vectorCopy(_Parameter, particleList);  //把上一次计算的结果vec保存到v
        
        LLG(_Parameter, particleList,1);  //K1
        vadd((double)0.5, _Parameter, particleList, 1);
        
        double ft0 = _Parameter[0].ft;
        LLG(_Parameter, particleList, 2);  //K2
        vadd((double)0.5, _Parameter, particleList, 2);

        LLG(_Parameter, particleList, 3);  //K3
        vadd(1, _Parameter, particleList, 3);

        LLG(_Parameter, particleList, 4);  //K4

        vadd4(_Parameter, particleList);
        vadd(1, _Parameter, particleList, 4);
        
        //Debug.Log("vec.x" + particleList[5].vec.x);

        _Parameter[0].ft = ft0;
       
        //把particleList合成
    }

}