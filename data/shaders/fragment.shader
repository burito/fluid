#version 410 core
uniform sampler2D diffuse;

layout (location = 0) in vec3 fragNormal;
layout (location = 1) in vec2 fragUV;
layout (location = 2) in vec4 fragColor;

layout (location = 0) out vec4 outColor;

void main()
{

	outColor = texture( diffuse, fragUV);

	float ld = dot(fragNormal, vec3(0,-1,0));

	outColor = fragColor;
//	outColor.xyz = outColor.xyz * 0.3 + fragColor.xyz * ld;
//	outColor = vec4(fragNormal, 1.);
//	outColor = vec4(vec3(ld), 1.);
//	outColor = vec4(1.);
}