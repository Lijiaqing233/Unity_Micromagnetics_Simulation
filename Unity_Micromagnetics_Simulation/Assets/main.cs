using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System;
using System.Runtime.InteropServices;

using System.IO;

 


public class main : SimBase {
    public UiBase myUI;
    double deltaTime = 0.0f;
    public int nx, ny, nz;
    public double p_hxc, p_hyc, p_hzc, dx, dy, dz, dt, alpha;
    int ThreadBlockSize = 64;
    int blockPerGrid;
    ComputeBuffer ParticleBuffer, argsBuffer, ParameterBuffer;
    ComputeBuffer gridParticleHashBuffer, gridParticleIndexBuffer;
    private uint[] _args;
    [Range(0, 1)]
    private double lerpt;
    // [Range(0.001f, 0.008f)]
    private double size = 1f;  //0.01
    [SerializeField]
    private Mesh Particle_Mesh;


    //设置compute shader
    [SerializeField]
    ComputeShader _computeShader; 


    [SerializeField]
    private Material _material;
    [SerializeField]
    private Material CPU_material;

    //UI
    bool Init;
    bool StartBtn = true;
    bool ResetBtn = false;
    bool StopBtn = false;
    public bool useGpuBtn = true;
    public TextAsset saveTXT;

    bool isDouble = false;  //选择GPGPU的数据类型是double，还是float


    //sim 
    int simTime = 0; //计算计算次数->为了进行误差分析

    bool usingRk4 = true;  //select for usingRk4
    double ft = 0;
    int jstep = 0;   //time step
    double time = 0;
    double ee = 1e-9;  //収束係数
    int ipr = 40;


    //private List<Particle> particleList = new List<Particle>(); //リアルタイムシミュレーションのシミュレーションデータを格納
    private Parameter[] _Parameter = new Parameter[1];          //
    private Parameterf[] _Parameterfloat = new Parameterf[1];

    Particle[] myParticle;
    Particlef[] myParticleGPGPU;

    // Start is called before the first frame update
    private void Start() {
        //UI
        myUI.FindUI();

        //Init
        
        int n = nx*ny*nz;
        int sizeParameter;
        int sizeParticle;

        myParticle = new Particle[n+1];
        myParticle.Initialize();
        myParticleGPGPU = new Particlef[n+1];
        myParticleGPGPU.Initialize();

        if (isDouble == true)
        {
            sizeParameter = System.Runtime.InteropServices.Marshal.SizeOf(typeof(Parameter));
            sizeParticle = System.Runtime.InteropServices.Marshal.SizeOf(typeof(Particle));
        }
        else 
        {
            sizeParameter = System.Runtime.InteropServices.Marshal.SizeOf(typeof(Parameterf));
            sizeParticle = System.Runtime.InteropServices.Marshal.SizeOf(typeof(Particlef));
        }


        ComputeBuffer computeBuffer = new ComputeBuffer(n + 1, sizeParticle);
        ParticleBuffer = computeBuffer;

        Init = true;
        init(_Parameter, myParticle, nx, ny, nz, 0, 0, 0, dx, dy, dz, dt, 1, Init);

        initFloat(_Parameterfloat, myParticleGPGPU, nx, ny, nz, 0, 0, 0, (float)dx, (float)dy, (float)dz, (float)dt, 1, Init);
 

        //Init  buffer
       

        Debug.Log(blockPerGrid);

        //InitSim();     //Init Sim   Giving an initial state
         Init = false;  //既に一回目初期化した

        jstep = 0;
        init(_Parameter, myParticle, nx, ny, nz, p_hxc, p_hyc, p_hzc, dx, dy, dz, dt, alpha, Init);
        initFloat(_Parameterfloat, myParticleGPGPU, nx, ny, nz, (float)p_hxc, (float)p_hyc, (float)p_hzc, (float)dx, (float)dy, (float)dz, (float)dt, (float)alpha, Init);

        ParameterBuffer = new ComputeBuffer(1, sizeParameter);  //176

        if (isDouble == true)
        {
            ParticleBuffer.SetData(myParticle);
            ParameterBuffer.SetData(_Parameter);
        }
        else 
        {
            ParticleBuffer.SetData(myParticleGPGPU);
            ParameterBuffer.SetData(_Parameterfloat);
        }        

        _args = new uint[5] { 0, 0, 0, 0, 0 };
        argsBuffer = new ComputeBuffer(1, _args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
         //Thread config
        blockPerGrid = (_Parameter[0].n + ThreadBlockSize - 1) / ThreadBlockSize;
    

    }
    // Update is called once per frame




    private void Awake() {
        Application.targetFrameRate = 200;  //Set frame rate


    }

       

    private void Update() {
        
         deltaTime += (Time.unscaledDeltaTime - deltaTime) * 0.1f;
        myUI.UpdateDisplayText(_Parameter);
        if (ResetBtn == true)
        {
            ResetBtn = false;
            Start();
        }
        
        if (StopBtn == false)
        {  //stop btn

            if (StartBtn == true)
            {
                //***************Computer shader base code*******************************//
                if (useGpuBtn == true)
                {
                    Rk4_CS();
                  
                    //rk4(_Parameter, myParticle);
                    /*
                    Particle[] receiveArr = new Particle[nx*ny*nz+1];
                    Particlef[] receiveArrf = new Particlef[nx*ny*nz+1];
                    //Debug.Log("GPU Fps:" + );
                    if (isDouble == true)
                    {
                        ParticleBuffer.GetData(receiveArr);
                    }
                    else 
                    {
                        ParticleBuffer.GetData(receiveArrf);
                    }
                     
                    ++jstep;
                    time = jstep * _Parameter[0].dt;
                  

                    if (jstep == jstep / ipr * ipr)
                       {
                        if ((jstep % ipr) == 0)
                        { 
                            string path = Application.persistentDataPath + "compute shadrERR.txt"; 
                            path = "C:/Users/85707/Desktop/compute shadrERR.txt";
                            if (isDouble == true)
                            {
                               // Debug.Log(check(myParticle, receiveArr));
                                double rus = (check(myParticle, receiveArr));
                                File.AppendAllText(path, rus.ToString()+"\r\n");
                                outputResults(receiveArr, ft, time);

                            }
                            else
                            { 
                                /*Debug.Log(checkFloat(myParticle, receiveArrf));
                                double rus = (checkFloat(myParticle, receiveArrf)); 
                                File.AppendAllText(path, rus.ToString()+"\r\n");
                      //          outputResultsFloat(receiveArrf, ft, time);*/
                      //      }
                    //    }
                    //   }
                    
                   
                   
                    



                    argsBuffer.SetData(_args);
                    Graphics.DrawMeshInstancedIndirect(Particle_Mesh, 0, _material, new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)), argsBuffer);
                }
                //***************Computer shader base code*******************************//

                //***********CPU Base code*************************//
                else
                {

                    rk4(_Parameter, myParticle);


                    //Debug.Log("time(ms):" + deltime);
                    //_Parameter[0].ft = Math.Sqrt(_Parameter[0].ft / (_Parameter[0].n)) * _Parameter[0].m;
                    outputResults(myParticle, _Parameter[0].ft, time);
                    rendering(myParticle, _Parameter[0].n);
                    if (_Parameter[0].ft < ee)
                    {
                        Debug.Log("収束");
                        //break;

                    }

                }
                //***********CPU Base code*************************//

            }
        }
        else
        {                       //Render even if the calculation stops
            if (useGpuBtn == true)
            {
                Graphics.DrawMeshInstancedIndirect(Particle_Mesh, 0, _material, new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)), argsBuffer);
            }
            else
            {
                rendering(myParticle, _Parameter[0].n);
            }
        }



 }

       

    

    private void rendering(Particle[] myParticle, int size)
    {
        Matrix4x4[] cpuMatrix = new Matrix4x4[size+1];
        for (int i = 1; i <= size; i++)
        {
            Quaternion tempQ = new Quaternion((float)myParticle[i].vec.x, (float)myParticle[i].vec.y, (float)myParticle[i].vec.z, 0);

            Vector3 tempVector = new Vector3((float)myParticle[i].position.x, (float)myParticle[i].position.y, (float)myParticle[i].position.z);
            int _Size = 1;
            double alpha = (double)Math.Acos(myParticle[i].vec.z);  //pitch  x  <-vec.z
            double beta = (double)Math.Atan(myParticle[i].vec.x + myParticle[i].vec.y);//yaw    y  <-vec.y
            double Gamma = 1;  //roll   z  <-vec.x

            //把模拟的结果转换为矩阵
            Matrix4x4 tempMatrix = new Matrix4x4()
            {  
                m00 = (float)(_Size * Math.Cos(Gamma) * Math.Cos(beta)),
                m01 = (float)(-Math.Sin(Gamma) * Math.Cos(alpha) + Math.Cos(Gamma) * Math.Sin(beta) * Math.Sin(alpha)),
                m02 = (float)(Math.Sin(Gamma) * Math.Sin(alpha)) + (float)(Math.Cos(Gamma) * Math.Sin(beta) * Math.Cos(alpha)),
                m03 = (float)myParticle[i].position.x,
                m10 = (float)(Math.Sin(Gamma) * Math.Cos(beta)),
                m11 = (float)(_Size * (Math.Cos(Gamma) * Math.Cos(alpha) + Math.Sin(Gamma) * Math.Sin(beta) * Math.Sin(alpha))),
                m12 = (float)(-Math.Cos(Gamma) * Math.Sin(alpha) + Math.Sin(Gamma) * Math.Sin(beta) * Math.Cos(alpha)),
                m13 = (float)myParticle[i].position.y,
                m20 = (float)(-Math.Sin(beta)),
                m21 = (float)(Math.Cos(beta) * Math.Sin(alpha)),
                m22 = (float)(_Size * Math.Cos(beta) * Math.Cos(alpha)),
                m23 = (float)myParticle[i].position.z,
                m30 = 0,
                m31 = 0,
                m32 = 0,
                m33 = 1
            };
            cpuMatrix[i] = tempMatrix;
           // Graphics.DrawMesh(Particle_Mesh, tempMatrix, CPU_material, 0);
        }

        Graphics.DrawMeshInstanced(Particle_Mesh, 0, CPU_material, cpuMatrix, size);

    }


    private void outputResults(Particle[] particleList, double ft, double time) {
        //double smx = 0, smy = 0, smz = 0;
        double3 myVecAvr;
        int ipr = 40;
        size = _Parameter[0].n;

        time = jstep * _Parameter[0].dt;
      
        myVecAvr = VecAvr(particleList);

      
                Debug.Log("    time = " + time * 1e9 + "   ns");
        //Debug.Log("    f:   " + ft + "    x:" + myVecAvr.x + "    y:" + myVecAvr.y + "    z:" + myVecAvr.z);
        Debug.Log("    x:" + myVecAvr.x + "    y:" + myVecAvr.y + "    z:" + myVecAvr.z);

    }


    private void outputResultsFloat(Particlef[] particleList, double ft, double time)
    {
        //double smx = 0, smy = 0, smz = 0;
        float3 myVecAvr;
        int ipr = 40;
        size = _Parameter[0].n;

        time = jstep * _Parameter[0].dt;

        myVecAvr = VecAvrFloat(particleList);


        Debug.Log("    time = " + time * 1e9 + "   ns");
        Debug.Log("    f:   " + ft + "    x:" + myVecAvr.x + "    y:" + myVecAvr.y + "    z:" + myVecAvr.z);


    }


    private double3 VecAvr(Particle[] particleList)
    {
        //输入一个电磁粒子数组，计算Vec的粒子平均值
        double3 VecAvr;
        double smx = 0, smy = 0, smz = 0;

        for (int i = 1; i <= size; i++)
        {
            smx += particleList[i].vec.x;
            smy += particleList[i].vec.y;
            smz += particleList[i].vec.z;
        }
        VecAvr.x = smx / size;
        VecAvr.y = smy / size;
        VecAvr.z = smz / size;

        return VecAvr;

    }

    private float3 VecAvrFloat(Particlef[] particleList)
    {
        //输入一个电磁粒子数组，计算Vec的粒子平均值
        float3 VecAvr;
        float smx = 0, smy = 0, smz = 0;

        for (int i = 1; i <= size; i++)
        {
            smx += particleList[i].vec.x;
            smy += particleList[i].vec.y;
            smz += particleList[i].vec.z;
        }
        VecAvr.x = (float)(smx / size);
        VecAvr.y = (float)(smy / size);
        VecAvr.z = (float)(smz / size);

        return VecAvr;

    }



    private double check(Particle[] myParticle, Particle[] receiveArr) {
        //输入CPU-BASE和Compute shader-base模拟结果,进行误差分析  输出相对误差
        double VecCpu_avr = 0;
        double3 temp_F_G;   //保存F-G的值

        double3 check_VecAvr;
        double3 check_VecAvr_GPU;
        double Ei = 0;
        int jstep = 0;  //time step

            ++jstep;


        VecCpu_avr = 0;
        Ei = 0;
                    for (int i = 1; i <= _Parameter[0].n; i++)
                    { //计算VECavr   全部粒子的平均Vec
                        VecCpu_avr += Math.Sqrt(myParticle[i].vec.x * myParticle[i].vec.x
                            + myParticle[i].vec.y * myParticle[i].vec.y
                            + myParticle[i].vec.z * myParticle[i].vec.z);
                    }
                    VecCpu_avr = VecCpu_avr / _Parameter[0].n;


                    for (int i = 1; i <= _Parameter[0].n; i++)
                    {
                            
                        //Fi - Gi
                        temp_F_G.x = (Math.Abs(myParticle[i].vec.x - receiveArr[i].vec.x));
                        temp_F_G.y = (Math.Abs(myParticle[i].vec.y - receiveArr[i].vec.y));
                        temp_F_G.z = (Math.Abs(myParticle[i].vec.z - receiveArr[i].vec.z));
                        //|Fi-Gi|/Favr
                        Ei += Math.Sqrt(temp_F_G.x* temp_F_G.x + temp_F_G.y* temp_F_G.y + temp_F_G.z* temp_F_G.z)/ VecCpu_avr;
                        Ei /= _Parameter[0].n;
                        temp_F_G.x = 0; temp_F_G.y = 0; temp_F_G.x = 0;
                    }

                    //计算平均Vec
                    check_VecAvr_GPU = VecAvr(receiveArr);
                    check_VecAvr = VecAvr(myParticle);
                    

                    //Debug.Log(jstep*_Parameter[0].dt*1e9);
                   return Ei;
                  
          
    }


    private double checkFloat(Particle[] myParticle, Particlef[] receiveArr)
    {
        //输入CPU-BASE和Compute shader-base模拟结果,进行误差分析  输出相对误差
        double VecCpu_avr = 0;
        double3 temp_F_G;   //保存F-G的值

        double3 check_VecAvr;
        double3 check_VecAvr_GPU;
        double Ei = 0;
        int jstep = 0;  //time step

        ++jstep;


        VecCpu_avr = 0;
        Ei = 0;
        for (int i = 1; i <= _Parameter[0].n; i++)
        { //计算VECavr   全部粒子的平均Vec
            VecCpu_avr += Math.Sqrt(myParticle[i].vec.x * myParticle[i].vec.x
                + myParticle[i].vec.y * myParticle[i].vec.y
                + myParticle[i].vec.z * myParticle[i].vec.z);
        }
        VecCpu_avr = VecCpu_avr / _Parameter[0].n;


        for (int i = 1; i <= _Parameter[0].n; i++)
        {

            //Fi - Gi
            temp_F_G.x = (Math.Abs(myParticle[i].vec.x - receiveArr[i].vec.x));
            temp_F_G.y = (Math.Abs(myParticle[i].vec.y - receiveArr[i].vec.y));
            temp_F_G.z = (Math.Abs(myParticle[i].vec.z - receiveArr[i].vec.z));
            //|Fi-Gi|/Favr
            Ei += Math.Sqrt(temp_F_G.x* temp_F_G.x + temp_F_G.y* temp_F_G.y + temp_F_G.z* temp_F_G.z)/ VecCpu_avr;
            Ei /= _Parameter[0].n;
            temp_F_G.x = 0; temp_F_G.y = 0; temp_F_G.x = 0;
        }

        //计算平均Vec
        //check_VecAvr_GPU = VecAvr(receiveArr);
        //check_VecAvr = VecAvr(myParticle);


        //Debug.Log(jstep*_Parameter[0].dt*1e9);
        return Ei;


    }


    void InitSim(){             //  Giving an initial state
                               
        if (useGpuBtn == true)
        {
            //***************Computer shader base code*******************************//
            for (int i = 0; i < 10000; i++)
            {

                Rk4_CS();

                argsBuffer.SetData(_args);
                Graphics.DrawMeshInstancedIndirect(Particle_Mesh, 0, _material, new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)), argsBuffer);
            }
            //***************Computer shader base code*******************************//
        }

        else
        {
            //***********CPU Base code*************************//
            _Parameter[0].ft = 10000;//Init 
            while (_Parameter[0].ft > ee)
            {
                jstep++;
                time = jstep * _Parameter[0].dt;
                rk4(_Parameter, myParticle);
               // _Parameter[0].ft = Math.Sqrt(_Parameter[0].ft / (_Parameter[0].n)) * _Parameter[0].m;

                outputResults(myParticle, _Parameter[0].ft, time);

                if (jstep > 10000)
                {
                    Debug.Log("収束不能!");
                    break;
                }

                if (_Parameter[0].ft < ee)
                {
                    Debug.Log("収束！");
                }
            }

            //***********CPU Base code*************************//
        }
    }

    private void OnDestroy(){
      // particleList.Clear();
       //ParticleBuffer.Dispose();
      // ParameterBuffer.Dispose();
       //argsBuffer.Dispose();
    }

    void OnGUI()
    {
        int w = Screen.width, h = Screen.height;
        GUIStyle style = new GUIStyle();     
        style.alignment = TextAnchor.UpperLeft;
        style.fontSize = h * 3 / 100;
        //new Color (0.0f, 0.0f, 0.5f, 1.0f);
        style.normal.textColor = Color.white;
        double msec = deltaTime * 1000.0f;
        double fps = 1.0f / deltaTime;
        string textFPS = string.Format("{0:0.0} ms   |   ({1:0.} fps)", msec, fps);
        string textNx = string.Format("nx: {0:0.}   |   ny: {1:0.}   |  nz:  {2:0.})", _Parameter[0].nx, _Parameter[0].ny, _Parameter[0].nz);
        string textHxc = string.Format("Hxc: {0:0.0} |  Hyc: {1:0.0}  |   Hzc:  {2:0.0})", _Parameter[0].hxc, _Parameter[0].hyc, _Parameter[0].hzc);
        string textother = string.Format("Alpha: {0:0.0} |   Dt:  {1:0.0}   |   n:   {2:0.0})", _Parameter[0].al, _Parameter[0].dt * 1e9, _Parameter[0].n);
        string textusingGPU = (useGpuBtn == true) ? ("using GPU") : ("using CPU");
        
        string[] textAll = { textFPS, textNx, textHxc, textother, textusingGPU};
     
        int dy = 0;
        for (int i = 0; i < textAll.Length; i++)
        {   
            dy += 40;
            Rect rect = new Rect(0, dy, w, h * 2 / 100);
            GUI.Label(rect, textAll[i], style);
        }


    }

    
    private void Rk4_CS(){

        //与cs的buffer链接
        //运行第1pass
        //把上一次计算的vec保存
        int kernelId1 = _computeShader.FindKernel("CSMain");
        _computeShader.SetBuffer(kernelId1, "_ParticleBuffer", ParticleBuffer);   //链接cs的Buffer
        _computeShader.SetBuffer(kernelId1, "_ParameterBuffer", ParameterBuffer);
        _computeShader.Dispatch(kernelId1, blockPerGrid, 1, 1); //blockPerGrid    运行compute shader

        //buffer相关的标记
        _args[0] = (uint)Particle_Mesh.GetIndexCount(0);  //36
        _args[1] = (uint)_Parameter[0].n;
        _args[2] = (uint)Particle_Mesh.GetIndexStart(0);  //0
        _args[3] = (uint)Particle_Mesh.GetBaseVertex(0);  //0

        _material.SetBuffer("_ParticleBuffer", ParticleBuffer);  //compute shaderの計算結果をシェーダーに送る(リンク)。
    
        _material.SetMatrix("_GameobjectMatrix", transform.localToWorldMatrix);
        _material.SetFloat("_Size", 1);
    }



    public void Stop(){
        Debug.Log("stop!");
        if (StopBtn == false){
            StopBtn = true;
        }

        else{
            StopBtn = false;
        }
    }

    public void Reset(){
        Debug.Log("stop!");
        //reset btn
        //Garbage
      
       //particleList.Clear();
       ParticleBuffer.Dispose();
       ParameterBuffer.Dispose();
       argsBuffer.Dispose();
        Start();
    }

    private void StartB(){
        StartBtn = true;
    }

    private void useGPU(){
        useGpuBtn = true;
    }

    private  void useCPU(){
        useGpuBtn = false;
    }

    public void UpdateNx(){
        string Nx = myUI.Nx.text;
        nx = int.Parse(Nx);
    }
  
    public void updateNy(){
        string Ny = myUI.Ny.text;
        ny = int.Parse(Ny);
    }
    public void updateNz(){
        string Nz = myUI.Nz.text;
        nz = int.Parse(Nz);
    }
    public void updateHxc(){
        string hxc = myUI.Hxc.text;
        p_hxc = double.Parse(hxc);
    }
    public void updateHyc(){
        string hyc = myUI.Hyc.text;
        p_hyc = double.Parse(hyc);
    }
    public void updateHzc(){
        string hzc = myUI.Hzc.text;
        p_hzc = double.Parse(hzc);
    }

    public void updateAl(){
        string Al = myUI.Al.text;
        alpha = double.Parse(Al);
    }

    public void updateDt(){
        string Dt = myUI.Dt.text;
        dt = double.Parse(Dt);
    }

    public void useGpu() {
        if (useGpuBtn == true)
        {
            useGpuBtn = false;
        }
        else if (useGpuBtn == false) {
            useGpuBtn = true;
        }

    }

}

