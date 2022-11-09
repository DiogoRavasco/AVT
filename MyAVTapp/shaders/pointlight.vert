#version 330

uniform mat4 m_pvm;
uniform mat4 m_viewModel;
uniform mat4 m_Model;   //por causa do cubo para a skybox
uniform mat3 m_normal;
uniform mat4 m_View;

uniform vec4 directional;

uniform vec4 point1;
uniform vec4 point2;
uniform vec4 point3;
uniform vec4 point4;
uniform vec4 point5;
uniform vec4 point6;

uniform vec4 left;
uniform vec4 right;

uniform int reflect_perFrag; //reflect vector calculated in the frag shader

uniform int texMode;

in vec4 position;
in vec4 normal;    //por causa do gerador de geometria
in vec4 texCoord;
in vec4 tangent;

out Data {
	vec3 normal;
	vec3 eye;
	vec3 lightDir[9];
	vec2 tex_coord;
	vec3 skyboxTexCoord;
	vec3 reflected;
} DataOut;

void main () {
	vec3 n, t, b;
	vec3 lightDir[9], eyeDir;
	vec3 aux;

	vec4 pos = m_viewModel * position;

	DataOut.skyboxTexCoord = vec3(m_Model * position);	//Transforma��o de modela��o do cubo unit�rio 
	DataOut.skyboxTexCoord.x = - DataOut.skyboxTexCoord.x; //Texturas mapeadas no interior logo negar a coordenada x

	DataOut.tex_coord = texCoord.st;
	
	if(texMode == 1)  {  //convert eye and light vectors to tangent space

		//Calculate components of TBN basis in eye space
		t = normalize(m_normal * tangent.xyz);  
		b = tangent.w * cross(n,t);

		lightDir[0] = vec3(directional);
		lightDir[1] = vec3(point1 - pos);
		lightDir[2] = vec3(point2 - pos);
		lightDir[3] = vec3(point3 - pos);
		lightDir[4] = vec3(point4 - pos);
		lightDir[5] = vec3(point5 - pos);
		lightDir[6] = vec3(point6 - pos);
	
		lightDir[7] = vec3(left - pos);
		lightDir[8] = vec3(right - pos);

		for (int i = 0; i < 9; i++) {
			aux.x = dot(lightDir[0], t);
			aux.y = dot(lightDir[0], b);
			aux.z = dot(lightDir[0], n);
			
			lightDir[0] = normalize(aux);
			
			DataOut.lightDir[i] = lightDir[i];
		}

		eyeDir = vec3(-pos);

		aux.x = dot(eyeDir, t);
		aux.y = dot(eyeDir, b);
		aux.z = dot(eyeDir, n);
		eyeDir = normalize(aux);

		DataOut.normal = n;
		DataOut.eye = eyeDir;

	}

	else{

		DataOut.normal = normalize(m_normal * normal.xyz);
		DataOut.eye = vec3(-pos);
		DataOut.lightDir[0] = vec3(directional);

	

		DataOut.lightDir[1] = vec3(point1 - pos);
		DataOut.lightDir[2] = vec3(point2 - pos);
		DataOut.lightDir[3] = vec3(point3 - pos);
		DataOut.lightDir[4] = vec3(point4 - pos);
		DataOut.lightDir[5] = vec3(point5 - pos);
		DataOut.lightDir[6] = vec3(point6 - pos);
	
		DataOut.lightDir[7] = vec3(left - pos);
		DataOut.lightDir[8] = vec3(right - pos);
	}

	if((texMode == 7) && (reflect_perFrag == 0)) {  //calculate here the reflected vector
		DataOut.reflected = vec3 (transpose(m_View) * vec4 (vec3(reflect(-DataOut.eye, DataOut.normal)), 0.0)); //reflection vector in world coord
		DataOut.reflected.x= -DataOut.reflected.x; // as texturas foram mapeadas no interior da skybox 
	}

	gl_Position = m_pvm * position;	
}