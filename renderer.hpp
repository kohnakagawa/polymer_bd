#pragma once
#include <GLFW/glfw3.h>
#include <GL/glu.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include "mvector3.hpp"
#include "parameter.hpp"

/*
  Wrapper class for rendering
 */
class PolymerRenderer{
  GLFWwindow* ptr_window = nullptr;

  GLUquadric* sphere = nullptr;
  const GLfloat color[3] = {0.000, 0.500, 0.000};
  const GLfloat rad = 0.5;
  const int div_phi = 10, div_the = 10;
  
  
  static void error_callback(int error, const char* description){
    fputs(description, stderr);
  }
  
  static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
      glfwSetWindowShouldClose(window, GL_TRUE);
  }

  inline void RenderSphere(const float3& pos){
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
    glPushMatrix();
    glTranslated(pos.x, pos.y, pos.z);
    gluSphere(sphere, rad, div_the, div_phi);
    glPopMatrix();
  }
  
public:
  PolymerRenderer(){
    glfwSetErrorCallback(error_callback);
    
    if(!glfwInit()){
      std::cerr << "glfwInit is failed \n";
      std::exit(EXIT_FAILURE);
    }

    sphere = gluNewQuadric();
    gluQuadricNormals(sphere, GL_SMOOTH);
  }

  ~PolymerRenderer(){
    //clean up
    glfwDestroyWindow(ptr_window);
    glfwTerminate();
  }
  
  void OpenWindow(const int hei, const int wid,
		  const std::string& title)
  {
    ptr_window = glfwCreateWindow(hei, wid, title.c_str(), nullptr, nullptr);
    if(!ptr_window){
      glfwTerminate();
      std::cerr << "Cannot create window \n";
      std::exit(EXIT_FAILURE);
    }
    
    glfwMakeContextCurrent(ptr_window);
  }

  void SetSwapInterval(const int i){
    glfwSwapInterval(i);
  }

  void SetKeyCallbackFunc(){
    glfwSetKeyCallback(ptr_window, key_callback);
  }
  
  void RenderPolymer(const double3* pos){
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      const float3 pos_flt(pos[i].x, pos[i].y, pos[i].z);
      RenderSphere(pos_flt);
    }
  }
  
  bool WindowShouldBeClosed(){ return glfwWindowShouldClose(ptr_window); }
  
  inline void GetFlameSize(int& wid, int& hei){ glfwGetFramebufferSize(ptr_window, &wid, &hei); }
  
  inline void SwapBuffer(){ glfwSwapBuffers(ptr_window); }
  
  inline void PollEvents(){ glfwPollEvents(); }

  void Resize() {
    int wid = -1, hei = -1;
    glfwGetFramebufferSize(ptr_window, &wid, &hei);
    
    glViewport(0, 0, wid, hei);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
  }
};
