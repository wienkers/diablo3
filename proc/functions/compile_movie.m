function compile_movie(f_in,f_out,dT_movie)

% Estimate framerate --> ~3s per period
  % dt in inertial units
framerate = max(3, ceil(2*pi/dT_movie / 3))

delete(['./movie_' f_out  '.mp4']);
system(['ffmpeg -r ' num2str(framerate) ' -f image2 -s 1280x720 -i ./' f_in '%04dlres.png -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" ./movie_' f_out  '.mp4'])


end
