# 🚀 BioScribe AI Deployment Guide

## Quick Start (Development)

### Prerequisites
- Node.js 18+ and npm
- Python 3.8+ and pip
- Git

### 1. Frontend Setup
```bash
cd bioscribe-ai
npm install
npm run dev
```
Frontend runs on: `http://localhost:3000`

### 2. Backend Setup
```bash
cd backend
pip install fastapi uvicorn biopython motor pymongo
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
```
Backend API runs on: `http://localhost:8000`
API Documentation: `http://localhost:8000/docs`

## Production Deployment

### Frontend (Vercel)
1. Push code to GitHub
2. Connect repository to Vercel
3. Set environment variables:
   ```
   NEXT_PUBLIC_API_URL=https://your-backend-url.com
   ```
4. Deploy automatically

### Backend (Railway/Render)
1. Create `Procfile`:
   ```
   web: uvicorn main:app --host 0.0.0.0 --port $PORT
   ```
2. Set environment variables:
   ```
   MONGODB_URL=mongodb://your-mongodb-url
   DATABASE_NAME=bioscribe_ai
   ```
3. Deploy with requirements.txt

### Docker Deployment

#### Frontend Dockerfile
```dockerfile
FROM node:18-alpine
WORKDIR /app
COPY package*.json ./
RUN npm ci --only=production
COPY . .
RUN npm run build
EXPOSE 3000
CMD ["npm", "start"]
```

#### Backend Dockerfile
```dockerfile
FROM python:3.9-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
EXPOSE 8000
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
```

#### Docker Compose
```yaml
version: '3.8'
services:
  frontend:
    build: ./bioscribe-ai
    ports:
      - "3000:3000"
    environment:
      - NEXT_PUBLIC_API_URL=http://backend:8000
    depends_on:
      - backend

  backend:
    build: ./backend
    ports:
      - "8000:8000"
    environment:
      - MONGODB_URL=mongodb://mongo:27017
      - DATABASE_NAME=bioscribe_ai
    depends_on:
      - mongo

  mongo:
    image: mongo:5
    ports:
      - "27017:27017"
    volumes:
      - mongo_data:/data/db

volumes:
  mongo_data:
```

## Environment Variables

### Frontend (.env.local)
```
NEXT_PUBLIC_API_URL=http://localhost:8000
```

### Backend (.env)
```
MONGODB_URL=mongodb://localhost:27017
DATABASE_NAME=bioscribe_ai
```

## Health Checks

### Frontend
- Check: `http://localhost:3000`
- Should show BioScribe AI interface

### Backend
- Health: `http://localhost:8000/health`
- API Docs: `http://localhost:8000/docs`
- Status: `http://localhost:8000/`

## Troubleshooting

### Common Issues

1. **Backend won't start**
   ```bash
   pip install fastapi uvicorn biopython motor pymongo
   ```

2. **Frontend build fails**
   ```bash
   rm -rf node_modules package-lock.json
   npm install
   ```

3. **API connection fails**
   - Check NEXT_PUBLIC_API_URL environment variable
   - Verify backend is running on correct port
   - Check CORS settings in backend

4. **3D visualization not working**
   - Ensure 3Dmol.js CDN is accessible
   - Check browser console for errors
   - Verify WebGL support

### Performance Optimization

1. **Frontend**
   - Enable Next.js image optimization
   - Use dynamic imports for heavy components
   - Implement proper caching strategies

2. **Backend**
   - Use async/await properly
   - Implement connection pooling for MongoDB
   - Add request rate limiting
   - Use background tasks for heavy computations

3. **Database**
   - Index frequently queried fields
   - Implement data cleanup for old sessions
   - Use MongoDB Atlas for production

## Monitoring

### Recommended Tools
- **Frontend**: Vercel Analytics, Sentry
- **Backend**: FastAPI metrics, Prometheus
- **Database**: MongoDB Compass, Atlas monitoring
- **Uptime**: UptimeRobot, Pingdom

### Key Metrics
- API response times
- Error rates
- User session duration
- Protein analysis completion rates
- Drug generation success rates

## Security

### Production Checklist
- [ ] Enable HTTPS
- [ ] Set up proper CORS policies
- [ ] Implement rate limiting
- [ ] Add input validation
- [ ] Use environment variables for secrets
- [ ] Enable database authentication
- [ ] Set up monitoring and alerting
- [ ] Regular security updates

## Scaling

### Horizontal Scaling
- Use load balancers for multiple backend instances
- Implement Redis for session storage
- Use CDN for static assets
- Consider microservices architecture

### Vertical Scaling
- Increase server resources
- Optimize database queries
- Use caching strategies
- Implement background job queues
