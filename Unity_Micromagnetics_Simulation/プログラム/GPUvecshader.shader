// Upgrade NOTE: replaced '_Object2World' with 'unity_ObjectToWorld'
// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

Shader "shader test/particle"
{
	Properties
	{

		 //_Specular("Specular", Color) = (1,1,1,1)

		 //_Diffuse("Diffuse",Color) = (1,1,1,1)

		 //_Gloss("Gloss",Range(1,255)) = 10

	}
		SubShader
	{
		Tags { "RenderType" = "Opaque" }
		LOD 200
		Pass
		{
			CGPROGRAM
			//#pragma target 3
			#pragma vertex vert
			#pragma fragment frag

			#include "Lighting.cginc"
			#include "UnityCG.cginc"
			struct Particle
			{
				 float qxx,qyy,qzz,qxy,qxz,qyz,li,lj,lk;
				 float3 vec,b,v,K1,K2,K3,K4,position;
				 
			};

		   StructuredBuffer<Particle> _ParticleBuffer;
			fixed4 _Color;
			float _lerp,_Size;
			float4x4 _GameobjectMatrix;
			sampler2D _MainTex;
			float4 _MainTex_ST;

			fixed4 _Diffuse;
			fixed4 _Specular;
			float _Gloss;
			

			struct appdata {
				float4 vertex:POSITION;
				float3 normal : NORMAL;

			};

			struct v2f {
				float4 pos:SV_POSITION;
				float3 worldNormal:NORMAL;
				float3 worldPos:TEXCOORD1;
				float3 color:COLOR;
				float2 uv : TEXCOORD0;
			};

			float4x4 GetModelToWorldMatrix(float3 pos, float3 vec)
			{
				float pi = 3.1415926;
				float alpha = acos(vec.z);  //pitch  x  <-vec.z
				float beta = abs(normalize(atan2(vec.y, vec.x))) - 1;
				float Gamma = 1;  //roll   z  <-vec.x
				float4x4 transformMatrix = float4x4(
					_Size * cos(Gamma) * cos(beta),-sin(Gamma) * cos(alpha) + cos(Gamma) * sin(beta) * sin(alpha) , sin(Gamma) * sin(alpha) + cos(Gamma) * sin(beta) * cos(alpha),pos.x,
					sin(Gamma) * cos(beta), _Size * (cos(Gamma) * cos(alpha) + sin(Gamma) * sin(beta) * sin(alpha)), -cos(Gamma) * sin(alpha) + sin(Gamma) * sin(beta) * cos(alpha), pos.y,
					-sin(beta), cos(beta) * sin(alpha),                 _Size * cos(beta) * cos(alpha),                   pos.z,
					0,                 0,                  0,                       1
				);
				return transformMatrix;
			}

			v2f vert(appdata VS_IN,uint instanceID :SV_INSTANCEID)
			{
				v2f VS_OUT;
				Particle particle = _ParticleBuffer[instanceID + 1];
				float4x4 WorldMatrix = GetModelToWorldMatrix(particle.position.xyz, particle.v.xyz);
				WorldMatrix = mul(_GameobjectMatrix,WorldMatrix);
				VS_IN.vertex = mul(WorldMatrix, VS_IN.vertex);
				VS_OUT.pos = mul(UNITY_MATRIX_VP,VS_IN.vertex);
				
				VS_OUT.worldNormal = normalize(mul(VS_IN.normal, (float3x3)unity_WorldToObject));
				VS_OUT.worldPos = mul(unity_ObjectToWorld, VS_IN.vertex).xyz;
				float3 worldLight =  normalize(_WorldSpaceLightPos0.xyz);
				fixed3 AColor =  fixed3(0.5 * particle.v.x +0.5 , 0.5*particle.v.y +0.5, 0.5 * particle.v.z + 0.5)*1;
				fixed3 diffuse =  AColor * saturate(dot(VS_OUT.worldNormal, worldLight));
				VS_OUT.color = ((AColor + diffuse));  //phong shading
				return VS_OUT;
			}
			fixed4 frag(v2f FS_IN) :SV_Target
			{
			   return fixed4(FS_IN.color, 0);
			}
			ENDCG
		}

	}
		FallBack Off
}