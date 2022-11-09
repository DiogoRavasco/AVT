#version 330

out vec4 colorOut;

uniform sampler2D texmap;
uniform sampler2D texmap1;
uniform sampler2D texmap2;
uniform sampler2D texmap3;
uniform sampler2D texmap4;
uniform samplerCube cubeMap;
uniform sampler2D normalMap;
uniform sampler2D texNormalMap;

uniform int texMode;

uniform mat4 m_View;
uniform int reflect_perFrag; //reflect vector calculated in the frag shader

struct Materials {
	vec4 diffuse;
	vec4 ambient;
	vec4 specular;
	vec4 emissive;
	float shininess;
	int texCount;
};

uniform Materials mat;

uniform vec4 lampDirection;

uniform bool fog;

uniform bool directionalOn;
uniform bool pointsOn;
uniform bool spotsOn;
uniform bool shadowMode;
uniform bool bumpmap1;

in Data {
	vec3 normal;
	vec3 eye;
	vec3 lightDir[9];
	vec2 tex_coord;
	vec3 skyboxTexCoord;
	vec3 reflected;
} DataIn;

const float density = 0.1;
const float gradient = 1.5;

const float linear = 0.1;
const float expo = 0.1;

const float spotLightApeture = 0.7;
const float spotExp = 80.0;

const float reflect_factor = 0.9;

const vec4 shadowC = vec4(.6, .6, .6, 1);

void main() {
	
	if(shadowMode) {
		colorOut = shadowC;
		return;
	} 
	
	colorOut = vec4(0);
	
	if(!directionalOn && !pointsOn && !spotsOn) return;

	
	vec4 texel, texel1, cube_texel; 

	vec3 n;

	if(texMode == 1)  // lookup normal from normal map, move from [0,1] to [-1, 1] range, normalize
		n = normalize(2.0 * texture(texNormalMap, DataIn.tex_coord).rgb - 1.0);
	else
		n = normalize(DataIn.normal);

	vec3 e = normalize(DataIn.eye);

	
	for (int i = 0; i < 9; i++) {		if (i == 0 && !directionalOn) continue;
		if (i >= 1 && i < 7 && !pointsOn) continue;
		if (i >= 7 && !spotsOn) break; 


		float attenuation = 3;
		
		if ( i != 0) {
			float distanceToLight = length(DataIn.lightDir[i]);
			attenuation += linear * distanceToLight + expo * distanceToLight * distanceToLight;	
		}

		vec3 l = normalize(DataIn.lightDir[i]);

		float intensity = max(dot(n,l), 0.0);

		if (i >= 7) {
			attenuation +=1;
			vec3 spotLightDirection = normalize(vec3(-lampDirection));

			float spotCos = dot(spotLightDirection, l);

			if(spotCos < spotLightApeture) intensity = 0.0;

			else attenuation = 1 / pow(spotCos, spotExp);
		}

		vec4 spec = vec4(0.0);
	
		if (intensity > 0.0) {
			vec3 h = normalize(l + e);
			float intSpec = max(dot(h,n), 0.0);
			spec = mat.specular * pow(intSpec, mat.shininess);
		}

		colorOut += max(intensity * mat.diffuse + spec, mat.ambient) / attenuation;

		if(texMode == 0) // modulate diffuse color with texel color
		{
			texel = texture(texmap, DataIn.tex_coord);  // texel from lighwood.tga
			colorOut += max(intensity * mat.diffuse * texel + spec,0.07 * texel) / attenuation;
		}
		if(texMode == 1) // bump map
		{
			texel = texture(texmap, DataIn.tex_coord);  // texel from stone.tga
			
			if (bumpmap1) {
				colorOut = max((intensity * texel + spec) / 9, 0.2 * texel / 9);
			} else {
				colorOut -= max(intensity * mat.diffuse + spec, mat.ambient) / attenuation;
				colorOut += max((intensity * texel + spec) / attenuation, 0.2 * texel / attenuation);
			}
		}
		else if (texMode == 2) // diffuse color is replaced by texel color, with specular area or ambient (0.1*texel)
		{
			texel = texture(texmap2, DataIn.tex_coord);  // texel from stone.tga
			colorOut += max(intensity*texel + spec, 0.07*texel) / attenuation;
		}
		else if (texMode == 3) 
		{
			texel = texture(texmap3, DataIn.tex_coord);  
		
			if(texel.a == 0.0) discard;
			else
				colorOut = 0.7*texel;
			
		}
		else if (texMode == 4) // Particles
		{
			texel = texture(texmap4, DataIn.tex_coord);
			colorOut = 0.7 * texel;
		}
		else if (texMode == 5) //SkyBox
		{
			colorOut = texture(cubeMap, DataIn.skyboxTexCoord);
		}
		else if(texMode == 7) // Environmental cube mapping
		{
			if(reflect_perFrag == 1) {  //reflected vector calculated here
				vec3 reflected1 = vec3 (transpose(m_View) * vec4 (vec3(reflect(-e, n)), 0.0)); //reflection vector in world coord
				reflected1.x= -reflected1.x;   
				cube_texel = texture(cubeMap, reflected1);
			}
			else
				cube_texel = texture(cubeMap, DataIn.reflected); //use interpolated reflected vector calculated in vertex shader

			texel = texture(texmap, DataIn.tex_coord);  // texel from lighwood.tga
			vec4 aux_color = mix(texel, cube_texel, reflect_factor);
			aux_color = max(intensity*aux_color + spec, 0.1*aux_color);
			//colorOut += vec4(aux_color.rgb, 1.0); 
		    colorOut = vec4(cube_texel.rgb, 1.0);
		}
		else if (texMode == 8) // Flare
		{
			texel = texture(texmap, DataIn.tex_coord);  //texel from element flare texture
			if((texel.a == 0.0)  || (mat.diffuse.a == 0.0) ) discard;
			else
				colorOut = mat.diffuse * texel;
		}
		else // multitexturing
		{
			texel = texture(texmap2, DataIn.tex_coord);  // texel from lighwood.tga
			texel1 = texture(texmap, DataIn.tex_coord);  // texel from stone.tga
			colorOut += max(intensity*texel*texel1 + spec, 0.07*texel*texel1) / attenuation;
		}
	}


	
	if (fog) {
		float distance = length(DataIn.eye) - 1;
		distance = max(0, distance);
		float visibility = exp(-pow(distance * density, gradient));
		visibility = clamp(visibility, 0.0, 1.0);

		colorOut = mix(vec4(0.5, 0.5, 0.5, 1), colorOut, visibility);
	}

	


}